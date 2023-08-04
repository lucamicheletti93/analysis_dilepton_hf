"""
Script for the training of ML models to be used in HF triggers
It requires hipe4ml and hipe4ml_converter to be installed:
  pip install hipe4ml
  pip install hipe4ml_converter
"""

import os
import sys
import argparse
import pickle
import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import yaml

from hipe4ml import plot_utils
from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler
from hipe4ml_converter.h4ml_converter import H4MLConverter

# pylint: disable=too-many-statements, too-many-branches, too-many-locals


def data_prep(config, ipt):
    """
    function for data preparation

    Parameters
    -----------------
    - config (dict): dictionary with config read from a yaml file
    - ipt (int): pt bin index
    """

    seed_split = config["data_prep"]["seed_split"]
    out_dir = config["output"]["directory"]
    suffix = config["output"]["out_labels"]["suffix"]
    leg_labels = [
        config["output"]["out_labels"]["bkg"],
        config["output"]["out_labels"]["prompt"],
        config["output"]["out_labels"]["nonprompt"]
    ]
    test_f = config["data_prep"]["test_fraction"]
    pt_min = config["data_prep"]["ptmins"][ipt]
    pt_max = config["data_prep"]["ptmaxs"][ipt]
    mass_leftsb = config["data_prep"]["bkg_sidebands"]["mass_left"]
    mass_rightsb = config["data_prep"]["bkg_sidebands"]["mass_right"]
    tree_name = config["data_prep"]["inputs"]["tree_name"]

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    hdl_prompt = TreeHandler(config["data_prep"]["inputs"]["prompt"], tree_name)
    hdl_nonprompt = TreeHandler(config["data_prep"]["inputs"]["nonprompt"], tree_name)
    hdl_bkg = TreeHandler(config["data_prep"]["inputs"]["bkg"], tree_name)

    pt_cut = f"{pt_min} < fPt < {pt_max}"
    mass_cut = f"({mass_leftsb[0]} < fM < {mass_leftsb[1]} or {mass_rightsb[0]} < fM < {mass_rightsb[1]})"

    hdl_prompt.apply_preselections(
        f"{pt_cut} and 0.99 < fOriginMcRec < 1.01", inplace=True)
    hdl_nonprompt.apply_preselections(
        f"{pt_cut} and 1.99 < fOriginMcRec < 2.01", inplace=True)
    hdl_bkg.apply_preselections(f"{pt_cut} and {mass_cut}", inplace=True)

    df_prompt = hdl_prompt.get_data_frame()
    df_nonprompt = hdl_nonprompt.get_data_frame()
    df_bkg = hdl_bkg.get_data_frame()

    n_prompt = len(df_prompt)
    n_nonprompt = len(df_nonprompt)
    n_bkg = len(df_bkg)
    print(f"\nNumber of available candidates in {pt_min} < pt < {pt_max}: \n     "
          f"prompt: {n_prompt}\n     nonprompt: {n_nonprompt}\n     bkg: {n_bkg}\n")

    n_cand_min = min([n_prompt, n_nonprompt, n_bkg])
    share = config["data_prep"]["class_balance"]["share"]
    if share == "equal":
        n_bkg = n_prompt = n_nonprompt = n_cand_min
    elif share == "all_signal":
        n_bkg = min(
            [n_bkg, (n_prompt + n_nonprompt) * config["data_prep"]["class_balance"]["bkg_factor"]])
    else:
        print(f"ERROR: class_balance option {share} not implemented")
        sys.exit()

    print(f"\nNumber of candidates used for training and test in {pt_min} < pt < {pt_max}: \n     "
          f"prompt: {n_prompt}\n     nonprompt: {n_nonprompt}\n     bkg: {n_bkg}\n")

    df_tot = pd.concat(
        [df_bkg[:n_bkg],
         df_prompt[:n_prompt],
         df_nonprompt[:n_nonprompt]],
        sort=True
    )

    labels_array = np.array([0]*n_bkg + [1]*n_prompt + [2]*n_nonprompt)
    if 0 < test_f < 1:
        train_set, test_set, y_train, y_test = train_test_split(
            df_tot, labels_array, test_size=test_f, random_state=seed_split
        )
    else:
        print("ERROR: test_fraction must belong to ]0,1[")
        sys.exit(0)

    train_test_data = [train_set, y_train, test_set, y_test]
    del df_tot

    df_list = [df_bkg, df_prompt, df_nonprompt]

    # _____________________________________________
    plot_utils.plot_distr(df_list, df_bkg.columns, 100, leg_labels,
                          figsize=(12, 7), alpha=0.3, log=True, grid=False, density=True)
    plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99,
                        top=0.96, hspace=0.55, wspace=0.55)
    plt.savefig(f"{out_dir}/DistributionsAll_{suffix}_pt{pt_min}_{pt_max}.pdf")
    plt.close("all")

    # _____________________________________________
    corr_matrix_fig = plot_utils.plot_corr(df_list, df_bkg.columns, leg_labels)
    for fig, lab in zip(corr_matrix_fig, leg_labels):
        plt.figure(fig.number)
        plt.subplots_adjust(left=0.2, bottom=0.25, right=0.95, top=0.9)
        fig.savefig(f"{out_dir}/CorrMatrix_{suffix}_{lab}_pt{pt_min}_{pt_max}.pdf")

    return train_test_data


def train(config, train_test_data, ipt):  # pylint: disable=too-many-locals
    """
    Function for the training

    Parameters
    -----------------
    - config: dictionary with config read from a yaml file
    - train_test_data: list with training and test data
    - ipt (int): pt bin index
    """

    out_dir = config["output"]["directory"]
    leg_labels = [
        config["output"]["out_labels"]["bkg"],
        config["output"]["out_labels"]["prompt"],
        config["output"]["out_labels"]["nonprompt"]
    ]
    suffix = config["output"]["out_labels"]["suffix"]
    pt_min = config["data_prep"]["ptmins"][ipt]
    pt_max = config["data_prep"]["ptmaxs"][ipt]
    n_classes = len(np.unique(train_test_data[3]))
    model_clf = xgb.XGBClassifier(use_label_encoder=False)
    training_vars = config["ml"]["training_vars"]
    hyper_pars = config["ml"]["hyper_pars"][ipt]
    model_hdl = ModelHandler(model_clf, training_vars, hyper_pars)

    # hyperparameters optimization
    if config["ml"]["hyper_pars_opt"]["activate"]:
        model_hdl.optimize_params_optuna(
            train_test_data,
            config["ml"]["hyper_pars_opt"]["hyper_par_ranges"],
            cross_val_scoring="roc_auc_ovo",
            timeout=config['ml']['hyper_pars_opt']['timeout'],
            n_jobs=config['ml']['hyper_pars_opt']['njobs'],
            n_trials=config['ml']['hyper_pars_opt']['ntrials'],
            direction='maximize'
        )
    else:
        model_hdl.set_model_params(hyper_pars)

    # train and test the model with the updated hyper-parameters
    y_pred_test = model_hdl.train_test_model(
        train_test_data,
        True,
        output_margin=config["ml"]["raw_output"],
        average=config["ml"]["roc_auc_average"],
        multi_class_opt=config["ml"]["roc_auc_approach"]
    )

    # save model
    pt_label = f"_pt{pt_min}_{pt_max}"
    if os.path.isfile(f"{out_dir}/ModelHandler_{suffix}{pt_label}.pickle"):
        os.remove(f"{out_dir}/ModelHandler_{suffix}{pt_label}.pickle")
    if os.path.isfile(f"{out_dir}/ModelHandler_onnx_{suffix}{pt_label}.onnx"):
        os.remove(f"{out_dir}/ModelHandler_onnx_{suffix}{pt_label}.onnx")

    model_hdl.dump_model_handler(f"{out_dir}/ModelHandler_{suffix}{pt_label}.pickle")
    model_conv = H4MLConverter(model_hdl)
    model_conv.convert_model_onnx(1)
    model_conv.dump_model_onnx(f"{out_dir}/ModelHandler_onnx_{suffix}{pt_label}.onnx")

    # plots
    # _____________________________________________
    plt.rcParams["figure.figsize"] = (10, 7)
    fig_ml_output = plot_utils.plot_output_train_test(
        model_hdl,
        train_test_data,
        80,
        config['ml']['raw_output'],
        leg_labels,
        True,
        density=True
    )

    if n_classes > 2:
        for fig, lab in zip(fig_ml_output, leg_labels):
            fig.savefig(f'{out_dir}/MLOutputDistr_{lab}_{suffix}{pt_label}.pdf')
    else:
        fig_ml_output.savefig(f'{out_dir}/MLOutputDistr_{suffix}{pt_label}.pdf')

    # _____________________________________________
    plt.rcParams["figure.figsize"] = (10, 9)
    fig_roc_curve = plot_utils.plot_roc(
        train_test_data[3],
        y_pred_test,
        None,
        leg_labels,
        config['ml']['roc_auc_average'],
        config['ml']['roc_auc_approach']
    )
    fig_roc_curve.savefig(f'{out_dir}/ROCCurveAll_{suffix}{pt_label}.pdf')
    pickle.dump(fig_roc_curve, open(
        f'{out_dir}/ROCCurveAll_{suffix}{pt_label}.pickle', 'wb'))

    # _____________________________________________
    plt.rcParams["figure.figsize"] = (12, 7)
    fig_feat_importance = plot_utils.plot_feature_imp(
        train_test_data[2][train_test_data[0].columns],
        train_test_data[3],
        model_hdl,
        leg_labels
    )
    n_plot = n_classes if n_classes > 2 else 1
    for i_fig, fig in enumerate(fig_feat_importance):
        if i_fig < n_plot:
            lab = leg_labels[i_fig] if n_classes > 2 else ''
            fig.savefig(f'{out_dir}/FeatureImportance_{lab}_{suffix}{pt_label}.pdf')
        else:
            fig.savefig(f'{out_dir}/FeatureImportanceAll_{suffix}{pt_label}.pdf')


def main(config):
    """
    Main function

    Parameters
    -----------------
    - config: dictionary with config read from a yaml file
    """
    for ipt, _ in enumerate(config["data_prep"]["ptmins"]):
        train_test_data = data_prep(config, ipt)
        if config["ml"]["dotraining"]:
            train(config, train_test_data, ipt)

    os._exit(0)  # pylint: disable=protected-access


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description="Arguments")
    PARSER.add_argument("config", metavar="text", default="config_training.yml",
                        help="config file for training")
    ARGS = PARSER.parse_args()

    with open(ARGS.config, "r") as yml_cfg:  # pylint: disable=bad-option-value
        CFG = yaml.load(yml_cfg, yaml.FullLoader)

    main(CFG)
