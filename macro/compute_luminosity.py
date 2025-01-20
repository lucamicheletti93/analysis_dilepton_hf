"""
Script for the preparation of efficiency maps and reflection templates
"""
import sys
import argparse
import ROOT
ROOT.gROOT.SetBatch(True)

def compute_luminosity(infile_names, labels, outfile_name, trigger):
    """
    Main function for luminosity calculation
    """

    if trigger not in ["fDiMuon", "fDiElectron"]:
        print(f"ERROR: wrong trigger mask {trigger}, please fix it!")
        sys.exit()

    if len(labels) != len(infile_names):
        print("ERROR: Please provide a label for each input file!")
        sys.exit()

    sigma_tvx = 59400000000 # pb
    lumi_tvx, lumi_tvx_afterbc = {}, {}
    hist_lumi_tvx = ROOT.TH1D("hist_lumi_tvx", ";;integrated luminosity (pb^{#minus1})",
                              len(labels) + 1, 0.5, len(labels) + 1.5)
    hist_lumi_tvx_afterbccuts = ROOT.TH1D("hist_lumi_tvx_afterbccuts",
                                          ";;integrated luminosity (pb^{#minus1})",
                                          len(labels) + 1, 0.5, len(labels) + 1.5)
    hist_lumi_tvx.GetXaxis().SetBinLabel(1, "all")
    hist_lumi_tvx_afterbccuts.GetXaxis().SetBinLabel(1, "all")
    for ilabel, label in enumerate(labels):
        hist_lumi_tvx.GetXaxis().SetBinLabel(2+ilabel, label)
        hist_lumi_tvx_afterbccuts.GetXaxis().SetBinLabel(2+ilabel, label)

    for ilabel, (label, infile_name) in enumerate(zip(labels, infile_names)):
        infile = ROOT.TFile.Open(infile_name)
        zorro_dir = infile.Get("hf-candidate-creator-2prong/Zorro")
        num_tvx_bcseltask = infile.Get("bc-selection-task/hCounterTVX")
        num_tvx_afterbc_bcseltask = infile.Get("bc-selection-task/hCounterTVXafterBCcuts")
        eff_bccuts = num_tvx_afterbc_bcseltask.Integral()/num_tvx_bcseltask.Integral()

        num_tvx = 0
        original_tvx, original_toi, selections = 0., 0., 0.
        for key in zorro_dir.GetListOfKeys():
            dir_run = zorro_dir.Get(key.GetName())
            original_tvx += dir_run.Get("InspectedTVX").GetBinContent(1)
            # For some reason AnalysedTriggersOfInterest has funny values, we use AnalysedTriggers
            # TODO: more investigations needed
            bin_toi = dir_run.Get("AnalysedTriggers").GetXaxis().FindBin(trigger)
            original_toi += dir_run.Get("AnalysedTriggers").GetBinContent(bin_toi)
            bin_selections = dir_run.Get("Selections").GetXaxis().FindBin(trigger)
            selections += dir_run.Get("Selections").GetBinContent(bin_selections)
        if selections > 0:
            num_tvx += original_tvx * original_toi / selections
        lumi_tvx[label] = num_tvx / sigma_tvx # pb-1
        lumi_tvx_afterbc[label] = lumi_tvx[label] * eff_bccuts
        hist_lumi_tvx.SetBinContent(ilabel+2, lumi_tvx[label])
        hist_lumi_tvx_afterbccuts.SetBinContent(ilabel+2, lumi_tvx_afterbc[label])

    hist_lumi_tvx.SetBinContent(1, sum(lumi_tvx.values()))
    hist_lumi_tvx_afterbccuts.SetBinContent(1, sum(lumi_tvx_afterbc.values()))

    outfile = ROOT.TFile(outfile_name, "recreate")
    hist_lumi_tvx.Write()
    hist_lumi_tvx_afterbccuts.Write()
    outfile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--infiles", "-i", nargs="+", type=str,
                        default=["AnalysisResults_22.root",
                                 "AnalysisResults_23.root",
                                 "AnalysisResults_24.root"],
                        help="input AnalysisResults files")
    parser.add_argument("--labels", "-l",  nargs="+", type=str,
                        default=["2022", "2023", "2024"],
                        help="labels corresponding to each input file")
    parser.add_argument("--outfile", "-o", metavar="text",
                        default="../data_shared/luminosity_2022_2023_2024_fDiMuon.root",
                        help="output file name")
    parser.add_argument("--trigger", "-t", metavar="text",
                        default="fDiMuon",
                        help="trigger mask, options: [fDiMuon, fDiElectron]")
    args = parser.parse_args()

    compute_luminosity(args.infiles, args.labels, args.outfile, args.trigger)
