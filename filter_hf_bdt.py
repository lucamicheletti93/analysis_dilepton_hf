"""
Simple script to filter derived AO2Ds for BDT trainings
"""

import argparse
import os
import pandas as pd
import yaml
import uproot

def filter_derived_ao2d(config):
    """
    Main function for filter

    Parameters
    ------------------
    - config (dict): dictionary with configs from yaml file
    """

    for input_type, input_names in config["inputs"].items():
        df_list = []
        if len(input_names) == 0:
            continue
        for input_name in input_names:
            with uproot.open(input_name) as infile:
                used_df = []
                for key in infile.keys():
                    if "O2hfcandd0lite" in key:
                        df_name = key.split(sep="/")[0]
                        if df_name not in used_df:
                            used_df.append(df_name)
                            df_list.append(infile[key].arrays(library="pd"))
        if len(df_list) > 1:
            output_df = pd.concat(df_list)
        else:
            output_df = df_list[0]
        # we remove the reflected signal and we split prompt/nonprompt
        d0_string = "(fCandidateSelFlag == 1 and fFlagMc >= 0)"
        d0bar_string = "(fCandidateSelFlag == 2 and fFlagMc <= 0)"
        origin = 0
        if input_type == "prompt":
            origin = 1
        elif input_type == "nonprompt":
            origin = 2
        output_df.query(f"({d0_string} or {d0bar_string}) and fOriginMcRec == {origin}",
                        inplace=True)
        selections = config['selections']
        if selections is not None:
            output_df.query(selections, inplace=True)
        output_df = output_df[config['columns_to_keep']]
        outdir = config['output']['directory']
        suffix = config['output']['suffix'][input_type]
        outfile = uproot.recreate(os.path.join(outdir, f"tree_d0_{input_type}{suffix}.root"))
        outfile["treeD0"] = output_df


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='Arguments')
    PARSER.add_argument('config', metavar='text', default='config_filter.yml',
                        help='input config yaml file name')
    ARGS = PARSER.parse_args()

    with open(ARGS.config, "r") as yml_cfg:  # pylint: disable=bad-option-value
        CFG = yaml.load(yml_cfg, yaml.FullLoader)

    filter_derived_ao2d(CFG)
