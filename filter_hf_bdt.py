"""
Simple script to filter derived AO2Ds for BDT trainings
"""

import argparse
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

    for input_type, input_name in config["inputs"].items():
        with uproot.open(input_name) as infile:
            df_list = []
            for key in infile.keys():
                if "O2hfcandd0lite" in key:
                    df_list.append(infile[key].arrays(library="pd"))
            output_df = pd.concat(df_list)
            outfile = uproot.recreate(f"tree_d0_{input_type}.root")
            outfile["treeD0"] = output_df


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='Arguments')
    PARSER.add_argument('config', metavar='text', default='config_filter.yml',
                        help='input config yaml file name')
    ARGS = PARSER.parse_args()

    with open(ARGS.config, "r") as yml_cfg:  # pylint: disable=bad-option-value
        CFG = yaml.load(yml_cfg, yaml.FullLoader)

    filter_derived_ao2d(CFG)
