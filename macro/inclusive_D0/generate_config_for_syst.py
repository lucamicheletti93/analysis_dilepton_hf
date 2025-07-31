"""
Simple script to generate_config_for_syst config files for systematic uncertainties
"""

import os
import argparse
import json
import yaml

def generate_configs(original_config, original_config_cutvar, n_loose, n_tight, loosest, tightest, ncls_its, ncls_tpc, outputdir, infile_trkcuts_data, infile_trkcuts_mc):
    """
    Main function for the generation of config files for systematics 
    """

    with open(original_config, "r") as yaml_file:
        cfg = yaml.safe_load(yaml_file)

    with open(original_config_cutvar, "r") as json_file:
        cfg_cutvar = json.load(json_file)

    orig_bdt_cuts = cfg["bdt_cuts"]["bkg"]
    orig_suffix = cfg["output"]["efficiencies"]["suffix"]
    orig_suffix_data = cfg["output"]["rawyields"]["suffix"]

    for ilooser in range(n_loose):
        outfile_name = ""
        if "yml" in original_config:
            outfile_name = original_config.replace(".yml", f"_loose{ilooser}.yml")
        elif "yaml" in original_config:
            outfile_name = original_config.replace(".yaml", f"_loose{ilooser}.yml")
        outfile_name = os.path.join(outputdir, outfile_name)

        bdt_cuts = []
        for ibdt, bdt_cut in enumerate(orig_bdt_cuts):
            bdt_cuts.append(bdt_cut + (ilooser + 1) * (loosest[ibdt] - bdt_cut) / n_loose)
        cfg["bdt_cuts"]["bkg"] = bdt_cuts
        cfg["bdt_cuts"]["nonprompt"] = []
        cfg["output"]["efficiencies"]["directory"] = "systematics"
        cfg["output"]["efficiencies"]["suffix"] = orig_suffix + f"_loose{ilooser}"
        cfg["output"]["rawyields"]["directory"] = "systematics"
        cfg["output"]["rawyields"]["suffix"] = orig_suffix_data + f"_loose{ilooser}"
        cfg_cutvar["output"]["directory"] = "systematics"
        cfg_cutvar["output"]["file"] = f"promptfrac{orig_suffix}_loose{ilooser}.root"

        with open(outfile_name, "w") as yaml_file:
            yaml.dump(cfg, yaml_file, default_flow_style=False)

        outfile_name_cutvar = original_config_cutvar.replace(".json", f"_loose{ilooser}.json")
        outfile_name_cutvar = os.path.join(outputdir, outfile_name_cutvar)
        cfg_cutvar["central_efficiency"]["inputdir"] = "systematics"
        cfg_cutvar["central_efficiency"]["inputfile"] = f"efficiencies_nocutnp{orig_suffix}_loose{ilooser}.root"
        with open(outfile_name_cutvar, "w") as json_file:
            json.dump(cfg_cutvar, json_file, indent=4)

    for itighter in range(n_tight):
        outfile_name = ""
        if "yml" in original_config:
            outfile_name = original_config.replace(".yml", f"_tight{itighter}.yml")
        elif "yaml" in original_config:
            outfile_name = original_config.replace(".yaml", f"_tight{itighter}.yml")
        outfile_name = os.path.join(outputdir, outfile_name)

        bdt_cuts = []
        for ibdt, bdt_cut in enumerate(orig_bdt_cuts):
            bdt_cuts.append(bdt_cut - (itighter + 1) * (bdt_cut - tightest[ibdt]) / n_tight)
        cfg["bdt_cuts"]["bkg"] = bdt_cuts
        cfg["bdt_cuts"]["nonprompt"] = []
        cfg["output"]["efficiencies"]["directory"] = "systematics"
        cfg["output"]["efficiencies"]["suffix"] = orig_suffix + f"_tight{itighter}"
        cfg["output"]["rawyields"]["directory"] = "systematics"
        cfg["output"]["rawyields"]["suffix"] = orig_suffix_data + f"_tight{itighter}"

        with open(outfile_name, "w") as yaml_file:
            yaml.dump(cfg, yaml_file, default_flow_style=False)

        outfile_name_cutvar = original_config_cutvar.replace(".json", f"_tight{itighter}.json")
        outfile_name_cutvar = os.path.join(outputdir, outfile_name_cutvar)
        cfg_cutvar["central_efficiency"]["inputdir"] = "systematics"
        cfg_cutvar["central_efficiency"]["inputfile"] = f"efficiencies_nocutnp{orig_suffix}_tight{itighter}.root"
        cfg_cutvar["output"]["directory"] = "systematics"
        cfg_cutvar["output"]["file"] = f"promptfrac{orig_suffix}_tight{itighter}.root"

        with open(outfile_name_cutvar, "w") as json_file:
            json.dump(cfg_cutvar, json_file, indent=4)

    if len(ncls_its) + len(ncls_tpc) > 0:
        itrkcut = 0
        for ncls_its_val in ncls_its:
            for ncls_tpc_val in ncls_tpc:
                
                outfile_name = ""
                if "yml" in original_config:
                    outfile_name = original_config.replace(".yml", f"_trkcut{itrkcut}.yml")
                elif "yaml" in original_config:
                    outfile_name = original_config.replace(".yaml", f"_trkcut{itrkcut}.yml")
                outfile_name = os.path.join(outputdir, outfile_name)

                cfg["trk_cuts"]["apply"] = True
                cfg["trk_cuts"]["ncls_its"] = ncls_its_val
                cfg["trk_cuts"]["ncls_tpc"] = ncls_tpc_val
                cfg["bdt_cuts"]["bkg"] = orig_bdt_cuts
                cfg["bdt_cuts"]["nonprompt"] = []
                cfg["output"]["efficiencies"]["directory"] = "systematics"
                cfg["output"]["efficiencies"]["suffix"] = orig_suffix + f"_trkcut{itrkcut}"
                cfg["output"]["rawyields"]["directory"] = "systematics"
                cfg["output"]["rawyields"]["suffix"] = orig_suffix_data + f"_trkcut{itrkcut}"

                cfg["input"]["data"] = infile_trkcuts_data
                cfg["input"]["mc"] = [infile_trkcuts_mc]
                cfg["input"]["corrbkg"] = [infile_trkcuts_mc]
                cfg["input"]["mc_weights"] = [1.]

                with open(outfile_name, "w") as yaml_file:
                    yaml.dump(cfg, yaml_file, default_flow_style=False)

                outfile_name_cutvar = original_config_cutvar.replace(".json", f"_trkcut{itrkcut}.json")
                outfile_name_cutvar = os.path.join(outputdir, outfile_name_cutvar)
                cfg_cutvar["central_efficiency"]["inputdir"] = "systematics"
                cfg_cutvar["central_efficiency"]["inputfile"] = f"efficiencies_nocutnp{orig_suffix}_trkcut{itrkcut}.root"
                cfg_cutvar["output"]["directory"] = "systematics"
                cfg_cutvar["output"]["file"] = f"promptfrac{orig_suffix}_trkcut{itrkcut}.root"

                with open(outfile_name_cutvar, "w") as json_file:
                    json.dump(cfg_cutvar, json_file, indent=4)

                itrkcut += 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--cfg_file", "-ce", metavar="text",
                        default="config_dzero_mb_y06.yml", help="original config file for raw yields and efficiency",
                        required=False)
    parser.add_argument("--config_file_cutvar", metavar="text",
                        default="config_cutvar_dzero_mb_y06.json", help="original config file for the cut-variation method",
                        required=False)
    parser.add_argument("--n_steps_loose", "-nl", type=int,
                        default=10, help="number of looser cut sets to be generate_eff_configsd",
                        required=False)
    parser.add_argument("--n_steps_tight", "-nt", type=int,
                        default=4, help="number of tighter cut sets to be generate_eff_configsd",
                        required=False)
    parser.add_argument("--ncls_its", "-nits", nargs="+", type=float,
                        default=[5., 6., 7.],
                        help="min number of ITS clusters to test", required=False)
    parser.add_argument("--ncls_tpc", "-ntpc", nargs="+", type=float,
                        default=[70., 80., 90., 100., 110.],
                        help="min number of TPC clusters to test", required=False)
    parser.add_argument("--loosest", "-l", nargs="+", type=float,
                        default=[0.15, 0.15, 0.15, 0.2, 0.3, 0.3, 0.4, 0.45, 0.6, 1.0, 1.0, 1.0, 1.0],
                        help="loosest cuts", required=False)
    parser.add_argument("--tightest", "-t", nargs="+", type=float,
                        default=[0.01, 0.01, 0.01, 0.02, 0.04, 0.04, 0.1, 0.15, 0.2, 0.1, 0.1, 0.1, 0.1],
                        help="tightest cuts", required=False)
    parser.add_argument("--infile_trkcuts_data", "-itd", metavar="text",
                        default="../../data_shared/2024_pass1_minbias/AnalysisResults_LHC24_minBias_sampled_D0_trksys.root",
                        help="input file for tracking systematics (data)", required=False)
    parser.add_argument("--infile_trkcuts_mc", "-itm", metavar="text",
                        default="../../data_shared/effmaps/AnalysisResults_LHC24k3_trackTuner_ptSmearing1p5_phiDep_trksys.root",
                        help="input file for tracking systematics (MC)", required=False)
    parser.add_argument("--outdir", "-d", metavar="text",
                        default="systematics", help="output directory", required=False)

    args = parser.parse_args()

    generate_configs(args.cfg_file, args.config_file_cutvar, args.n_steps_loose, args.n_steps_tight,
                     args.loosest, args.tightest, args.ncls_its, args.ncls_tpc, args.outdir, args.infile_trkcuts_data, args.infile_trkcuts_mc)

