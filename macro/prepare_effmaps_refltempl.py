"""
Script for the preparation of efficiency maps and reflection templates
"""
import argparse
import math
import yaml
import numpy as np
import ROOT
ROOT.gROOT.SetBatch(True)

def get_reflections(config): #pylint:disable=too-many-locals,too-many-statements,too-many-branches
    """
    Main function for the preparation of the reflection templates
    """

    ROOT.gStyle.SetPalette(ROOT.kRainBow)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadLeftMargin(0.18)
    ROOT.gStyle.SetPadRightMargin(0.16)
    ROOT.gStyle.SetTitleOffset(1.4, "z")

    with open(config, 'r') as yml:
        cfg = yaml.load(yml, yaml.FullLoader)
    pt_mins = cfg['pt_mins']
    pt_maxs = cfg['pt_maxs']
    bdtbkg_cuts = cfg['bdtbkg_cuts']
    infile_names = cfg["inputs"]["files"]
    input_weights = cfg["inputs"]["weights"]

    sparse_rec, sparse_gen = [], []
    n_events = []
    for ifile, infile_name in enumerate(infile_names):
        infile = ROOT.TFile.Open(infile_name)
        hist_events = infile.Get("hf-candidate-creator-2prong/hCollisions")
        n_events.append(hist_events.GetBinContent(1)) # to reweight
        sparse_rec.append(infile.Get("hf-task-d0/hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"))
        sparse_gen.append(infile.Get("hf-task-d0/hSparseAcc"))

    for ifile, nev in enumerate(n_events):
        weights_ev_mc = sum(n_events)/nev # we normalise all the inputs to the same weight
        input_weights[ifile] *= weights_ev_mc # we weight them according to the values in the config

    hist_refl, hist_signal = None, None
    hist_reco_p_toreb, hist_reco_np_toreb = None, None
    for ipt, (pt_min, pt_max, bdtbkg_cut) in enumerate(zip(pt_mins, pt_maxs, bdtbkg_cuts)):
        pt_bin_min = sparse_rec[0].GetAxis(4).FindBin(pt_min * 1.0001)
        pt_bin_max = sparse_rec[0].GetAxis(4).FindBin(pt_max * 0.9999)
        bdtbkg_bin = sparse_rec[0].GetAxis(0).FindBin(bdtbkg_cut * 0.9999)
        for ifile, _ in enumerate(infile_names):
            sparse_rec[ifile].GetAxis(4).SetRange(pt_bin_min, pt_bin_max)
            sparse_rec[ifile].GetAxis(0).SetRange(1, bdtbkg_bin)
            if ipt == 0 and ifile == 0:
                sparse_rec[ifile].GetAxis(6).SetRange(3, 4)
                hist_refl = sparse_rec[ifile].Projection(3)
                hist_refl.SetName("hist_refl")
                sparse_rec[ifile].GetAxis(6).SetRange(1, 2)
                hist_signal = sparse_rec[ifile].Projection(3)
                hist_signal.Scale(input_weights[ifile])
                hist_signal.SetName("hist_signal")
                # prompt reco pt vs y map
                sparse_rec[ifile].GetAxis(8).SetRange(2, 2)
                hist_reco_p_toreb = sparse_rec[ifile].Projection(5, 4, "O")
                hist_reco_p_toreb.Sumw2()
                hist_reco_p_toreb.Scale(input_weights[ifile])
                # non-prompt reco pt vs y map
                sparse_rec[ifile].GetAxis(8).SetRange(3, 3)
                hist_reco_np_toreb = sparse_rec[ifile].Projection(5, 4, "O")
                hist_reco_np_toreb.Sumw2()
                hist_reco_np_toreb.Scale(input_weights[ifile])
                sparse_rec[ifile].GetAxis(8).SetRange(-1, -1)
            else:
                sparse_rec[ifile].GetAxis(6).SetRange(3, 4)
                hist_refl.Add(sparse_rec[ifile].Projection(3))
                sparse_rec[ifile].GetAxis(6).SetRange(1, 2)
                hist_signal.Add(sparse_rec[ifile].Projection(3), input_weights[ifile])
                # prompt reco pt vs y map
                sparse_rec[ifile].GetAxis(8).SetRange(2, 2)
                hist_reco_p_toreb.Add(sparse_rec[ifile].Projection(5, 4, "O"), input_weights[ifile])
                # non-prompt reco pt vs y map
                sparse_rec[ifile].GetAxis(8).SetRange(3, 3)
                hist_reco_np_toreb.Add(sparse_rec[ifile].Projection(5, 4, "O"), input_weights[ifile])
                sparse_rec[ifile].GetAxis(8).SetRange(-1, -1)

    hist_gen_p_toreb, hist_gen_np_toreb = None, None
    for ifile, _ in enumerate(infile_names):
        if ifile == 0:
            # prompt gen pt vs y map
            sparse_gen[ifile].GetAxis(3).SetRange(2, 2)
            hist_gen_p_toreb = sparse_gen[ifile].Projection(2, 0, "O")
            hist_gen_p_toreb.Sumw2()
            hist_gen_p_toreb.Scale(input_weights[ifile])
            # non-prompt gen pt vs y map
            sparse_gen[ifile].GetAxis(3).SetRange(3, 3)
            hist_gen_np_toreb = sparse_gen[ifile].Projection(2, 0, "O")
            hist_gen_np_toreb.Sumw2()
            hist_gen_np_toreb.Scale(input_weights[ifile])
            sparse_gen[ifile].GetAxis(3).SetRange(-1, -1)
        else:
            # prompt gen pt vs y map
            sparse_gen[ifile].GetAxis(3).SetRange(2, 2)
            hist_gen_p_toreb.Add(sparse_gen[ifile].Projection(2, 0, "O"), input_weights[ifile])
            # non-prompt gen pt vs y map
            sparse_gen[ifile].GetAxis(3).SetRange(3, 3)
            hist_gen_np_toreb.Add(sparse_gen[ifile].Projection(2, 0, "O"), input_weights[ifile])
            sparse_gen[ifile].GetAxis(3).SetRange(-1, -1)

    # efficiency maps
    # rebin by hand rec to have compatible bins wrt gen, pT-variable rebin not supported for TH2
    pt_lims = [0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0,
               8.5, 9.0, 9.5, 10., 11., 12., 13., 14., 15., 16., 18., 20., 22., 24., 30., 36., 50.]
    nptbins = len(pt_lims) - 1
    nybins = hist_reco_p_toreb.GetYaxis().GetNbins()
    y_lims = [hist_reco_p_toreb.GetYaxis().GetBinLowEdge(iy) for iy in range(1, nybins + 1)]
    y_lims.append(hist_reco_p_toreb.GetYaxis().GetBinUpEdge(nybins))
    pt_lims = np.array(pt_lims, dtype=np.float64)
    y_lims = np.array(y_lims, dtype=np.float64)
    hist_reco_p = ROOT.TH2D("hist_recoD0_prompt", ";#it{p}_{T} (GeV/#it{c});#it{y}",
                            nptbins, pt_lims, nybins, y_lims)
    hist_reco_np = ROOT.TH2D("hist_recoD0_nonprompt", ";#it{p}_{T} (GeV/#it{c});#it{y}",
                             nptbins, pt_lims, nybins, y_lims)
    for iy in range(1, hist_reco_p.GetYaxis().GetNbins() + 1):
        for ipt in range(1, hist_reco_p.GetXaxis().GetNbins() + 1):
            ptmin = hist_reco_p.GetXaxis().GetBinLowEdge(ipt)
            ptmax = hist_reco_p.GetXaxis().GetBinUpEdge(ipt)
            counts_p, error_p, counts_np, error_np = 0., 0., 0., 0.
            for iptorig in range(1, hist_reco_p_toreb.GetXaxis().GetNbins() + 1):
                ptmin_orig = hist_reco_p_toreb.GetXaxis().GetBinLowEdge(iptorig)
                ptmax_orig = hist_reco_p_toreb.GetXaxis().GetBinUpEdge(iptorig)
                if (ptmin_orig >= ptmin or math.isclose(ptmin_orig, ptmin)) and \
                    (ptmax_orig <= ptmax or math.isclose(ptmax_orig, ptmax)):
                    counts_p += hist_reco_p_toreb.GetBinContent(iptorig, iy)
                    error_p += hist_reco_p_toreb.GetBinError(iptorig, iy)**2
                    counts_np += hist_reco_np_toreb.GetBinContent(iptorig, iy)
                    error_np += hist_reco_np_toreb.GetBinError(iptorig, iy)**2
            hist_reco_p.SetBinContent(ipt, iy, counts_p)
            hist_reco_p.SetBinError(ipt, iy, np.sqrt(error_p))
            hist_reco_np.SetBinContent(ipt, iy, counts_np)
            hist_reco_np.SetBinError(ipt, iy, np.sqrt(error_np))

    # rebin by hand gen (different pT bins), pT-variable rebin not supported for TH2
    hist_gen_p = hist_reco_p.Clone("hist_genD0_prompt")
    hist_gen_p.Reset()
    hist_gen_np = hist_reco_np.Clone("hist_genD0_nonprompt")
    hist_gen_np.Reset()
    for iy in range(1, hist_gen_p.GetYaxis().GetNbins() + 1):
        for ipt in range(1, hist_gen_p.GetXaxis().GetNbins() + 1):
            ptmin = hist_gen_p.GetXaxis().GetBinLowEdge(ipt)
            ptmax = hist_gen_p.GetXaxis().GetBinUpEdge(ipt)
            counts_p, error_p, counts_np, error_np = 0., 0., 0., 0.
            for iptorig in range(1, hist_gen_p_toreb.GetXaxis().GetNbins() + 1):
                ptmin_orig = hist_gen_p_toreb.GetXaxis().GetBinLowEdge(iptorig)
                ptmax_orig = hist_gen_p_toreb.GetXaxis().GetBinUpEdge(iptorig)
                if (ptmin_orig >= ptmin or math.isclose(ptmin_orig, ptmin)) and \
                    (ptmax_orig <= ptmax or math.isclose(ptmax_orig, ptmax)):
                    counts_p += hist_gen_p_toreb.GetBinContent(iptorig, iy)
                    error_p += hist_gen_p_toreb.GetBinError(iptorig, iy)**2
                    counts_np += hist_gen_np_toreb.GetBinContent(iptorig, iy)
                    error_np += hist_gen_np_toreb.GetBinError(iptorig, iy)**2
            hist_gen_p.SetBinContent(ipt, iy, counts_p)
            hist_gen_p.SetBinError(ipt, iy, np.sqrt(error_p))
            hist_gen_np.SetBinContent(ipt, iy, counts_np)
            hist_gen_np.SetBinError(ipt, iy, np.sqrt(error_np))

    hist_eff_p = hist_reco_p.Clone("histAxeD0_prompt")
    hist_eff_p.Sumw2()
    hist_eff_p.Divide(hist_reco_p, hist_gen_p, 1., 1., "B")
    hist_eff_np = hist_reco_np.Clone("histAxeD0_nonprompt")
    hist_eff_np.Sumw2()
    hist_eff_np.Divide(hist_reco_np, hist_gen_np, 1., 1., "B")

    canv = ROOT.TCanvas("canv", "", 1000, 500)
    canv.Divide(2, 1)
    canv.cd(1)
    hist_eff_p.GetZaxis().SetTitle("(Acc #times #font[152]{e}) prompt D^{0}")
    hist_eff_p.DrawCopy("colz")
    canv.cd(2)
    hist_eff_np.GetZaxis().SetTitle("(Acc #times #font[152]{e}) non-prompt D^{0}")
    hist_eff_np.DrawCopy("colz")

    outfile_name_eff = cfg["outputs"]["efficiencies"]
    canv.SaveAs(outfile_name_eff.replace(".root", ".pdf"))
    outfile = ROOT.TFile(outfile_name_eff, "recreate")
    canv.Write()
    hist_reco_p.Write()
    hist_reco_np.Write()
    hist_gen_p.Write()
    hist_gen_np.Write()
    hist_eff_p.Write()
    hist_eff_np.Write()
    outfile.Close()

    # reflections
    hist_r = ROOT.TH1F("hist_r", ";;R/S", 1, -0.5, 0.5)
    hist_s = ROOT.TH1F("hist_s", ";;R/S", 1, -0.5, 0.5)
    hist_r.Sumw2()
    hist_s.Sumw2()
    hist_r.SetBinContent(1, hist_refl.Integral())
    hist_s.SetBinContent(1, hist_signal.Integral())
    hist_r_over_s = hist_r.Clone("hist_r_over_s")
    hist_r_over_s.Divide(hist_s)

    hist_refl_smooth = hist_refl.Clone("hist_refl_smooth")
    hist_refl_smooth.Smooth(100)

    outfile_name_refl = cfg["outputs"]["reflections"]
    outfile = ROOT.TFile(outfile_name_refl, "recreate")
    hist_refl.Write()
    hist_refl_smooth.Write()
    hist_signal.Write()
    hist_r_over_s.Write()
    outfile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--config", "-c", metavar="text",
                        default="config_eff_and_refl.yml",
                        help="config file with BDT cuts")
    args = parser.parse_args()

    get_reflections(args.config)
