"""
Get integrated luminosity from AnalysisResults.root files
"""

import argparse
import ROOT


def get_lumi(infile_names, suffix):
    """

    """

    hist_lumi_per_run = []
    lumi = 0
    for infile_name in infile_names:
        infile = ROOT.TFile.Open(infile_name)
        hist_collisions = infile.Get("hf-candidate-creator-dstar/hCollisions")
        if not isinstance(hist_collisions, ROOT.TH1):
            hist_collisions = infile.Get("hf-candidate-creator-3prong/hCollisions")
        if not isinstance(hist_collisions, ROOT.TH1):
            hist_collisions = infile.Get("hf-candidate-creator-2prong/hCollisions")
        ibin_zvtx = hist_collisions.GetXaxis().FindBin("PV #it{z}")
        zvtx_eff = hist_collisions.GetBinContent(ibin_zvtx) / hist_collisions.GetBinContent(ibin_zvtx-1)
        hist_lumi_per_run.append(infile.Get("bc-selection-task/hLumiTVXafterBCcuts"))
        print(hist_lumi_per_run[-1].GetNbinsX())
        for ibin in range(1, hist_lumi_per_run[-1].GetNbinsX()+1):
            lumi += hist_lumi_per_run[-1].GetBinContent(ibin) * zvtx_eff # mub-1

    hist_lumi = ROOT.TH1F("hist_lumi", ";;luminosity (pb^{#minus1})", 1, 0., 1.)
    hist_lumi.SetBinContent(1, lumi)

    outfile = ROOT.TFile(f"luminosity{suffix}.root", "recreate")
    hist_lumi.Write()
    outfile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--infiles", "-i", nargs="+", type=str,
                        default=[], help="input AnalysisResults files")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output file")
    args = parser.parse_args()

    get_lumi(args.infiles, args.suffix)
