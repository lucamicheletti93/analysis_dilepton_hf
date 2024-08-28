"""
"""
import argparse
import yaml
import ROOT
ROOT.gROOT.SetBatch(True)

def get_reflections(infile_name, outfile_name, config):
    """
    Main function for the preparation of the reflection templates
    """

    with open(config, 'r') as yml:
        config_cuts = yaml.load(yml, yaml.FullLoader)
    pt_mins = config_cuts['pt_mins']
    pt_maxs = config_cuts['pt_maxs']
    bdtbkg_cuts = config_cuts['bdtbkg_cuts']

    infile = ROOT.TFile.Open(infile_name)
    sparse = infile.Get("hf-task-d0/hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type")
    sparse.GetAxis(2).SetRange(41, 240) # remove the large mass region, where we only have high-pT

    hist_refl, hist_signal = None, None
    for ipt, (pt_min, pt_max, bdtbkg_cut) in enumerate(zip(pt_mins, pt_maxs, bdtbkg_cuts)):
        pt_bin_min = sparse.GetAxis(3).FindBin(pt_min * 1.0001)
        pt_bin_max = sparse.GetAxis(3).FindBin(pt_max * 0.9999)
        bdtbkg_bin = sparse.GetAxis(0).FindBin(bdtbkg_cut * 0.9999)
        sparse.GetAxis(3).SetRange(pt_bin_min, pt_bin_max)
        sparse.GetAxis(0).SetRange(1, bdtbkg_bin)
        if ipt == 0:
            sparse.GetAxis(7).SetRange(3, 4)
            hist_refl = sparse.Projection(2)
            hist_refl.SetName("hist_refl")
            sparse.GetAxis(7).SetRange(1, 2)
            hist_signal = sparse.Projection(2)
            hist_signal.SetName("hist_signal")
        else:
            sparse.GetAxis(7).SetRange(3, 4)
            hist_refl.Add(sparse.Projection(2))
            sparse.GetAxis(7).SetRange(1, 2)
            hist_signal.Add(sparse.Projection(2))

    hist_r = ROOT.TH1F("hist_r", ";;R/S", 1, -0.5, 0.5)
    hist_s = ROOT.TH1F("hist_s", ";;R/S", 1, -0.5, 0.5)
    hist_r.SetBinContent(1, hist_refl.Integral())
    hist_s.SetBinContent(1, hist_signal.Integral())
    hist_r.Sumw2()
    hist_s.Sumw2()
    hist_r_over_s = hist_r.Clone("hist_r_over_s")
    hist_r_over_s.Divide(hist_s)

    hist_refl_smooth = hist_refl.Clone("hist_refl_smooth")
    hist_refl_smooth.Smooth(100)

    outfile = ROOT.TFile(outfile_name, "recreate")
    hist_refl.Write()
    hist_refl_smooth.Write()
    hist_signal.Write()
    hist_r_over_s.Write()
    outfile.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--input_file", "-i", metavar="text",
                        default="AnalysisResults.root",
                        help="input ROOT file")
    parser.add_argument("--output_file", "-o", metavar="text",
                        default="reflections_d0.root",
                        help="output ROOT file")
    parser.add_argument("--config", "-c", metavar="text",
                        default="config_bdt_cuts.yml",
                        help="config file with BDT cuts")
    args = parser.parse_args()

    get_reflections(args.input_file, args.output_file, args.config)
