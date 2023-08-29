"""
Script for the application of ML models for HF candidates
It requires hipe4ml and hipe4ml_converter to be installed:
  pip install hipe4ml
  pip install hipe4ml_converter
"""

import argparse
import yaml

from hipe4ml.tree_handler import TreeHandler
from hipe4ml.model_handler import ModelHandler

def apply(config):
    """
    Function for the training

    Parameters
    -----------------
    - config: dictionary with config read from a yaml file
    - use_onnx: use also onnx models
    """

    tree_handler = TreeHandler(config["inputs"]["data"], tree_name="treeD0")
    for ipt, (pt_min, pt_max) in enumerate(zip(config["ptmins"], config["ptmaxs"])):
        tree_handler_pt = tree_handler.apply_preselections(
            f"{pt_min} < fPt < {pt_max}", inplace=False)
        model_hdl = ModelHandler()
        model_hdl.load_model_handler(config["models"][ipt])
        tree_handler_pt.apply_model_handler(model_hdl, output_margin=False)
        tree_handler_pt.write_df_to_root_files(
            f"treeD0_data_applied_pt{pt_min:0.1f}_{pt_max:0.1f}", "treeD0")

if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description="Arguments")
    PARSER.add_argument("config", metavar="text", default="config_training.yml",
                        help="config file for training")
    ARGS = PARSER.parse_args()

    with open(ARGS.config, "r") as yml_cfg:  # pylint: disable=bad-option-value
        CFG = yaml.load(yml_cfg, yaml.FullLoader)

    apply(CFG)
