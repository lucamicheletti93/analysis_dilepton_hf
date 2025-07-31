"""
"""

import argparse
import uproot
import numpy as np
import ROOT

COLORS = [ROOT.kRed+1, ROOT.kAzure+4, ROOT.kGray+2]

def set_style(hist, color, markerstyle=ROOT.kFullCircle, marksersize=1):
    """
    """

    hist.SetLineWidth(2)
    hist.SetMarkerStyle(markerstyle)
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetMarkerSize(marksersize)


def produce_predictions(infile_name, outfile_name, kine_cuts, plot_lhcb):
    """
    """

    ROOT.gStyle.SetTitleSize(0.045, "xy")
    ROOT.gStyle.SetLabelSize(0.045, "xy")
    ROOT.gStyle.SetPadRightMargin(0.035)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadTopMargin(0.075)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

    infile_root = ROOT.TFile(infile_name)
    hist_xsec = infile_root.Get("histXsec_42")
    hist_numev = infile_root.Get("histNumEvents")
    xsec = hist_xsec.GetBinContent(1) * 1.e3 # mub
    nevents = hist_numev.GetBinContent(2)
    hist_jpsi = infile_root.Get("histPtVsY_443")
    hist_dmes = infile_root.Get("histPtVsY_421")
    hist_jpsi.SetDirectory(0)
    hist_dmes.SetDirectory(0)
    hist_dmes.Scale(1./2) # stored D0 + D0bar
    infile_root.Close()

    df = uproot.open(infile_name)["tuplePairs"].arrays(library="pd")

    cuts = ["no_kine_cuts"] + kine_cuts
    n_cuts = len(cuts)
    outfile = ROOT.TFile(outfile_name, "recreate")
    canv, leg = {}, {}
    hist_ratio_tot_jpsiincl_dmesincl = ROOT.TH1D("hist_ratio_tot_jpsiincl_dmesincl", ";;#it{R}_{DJ/#Psi} (#mub)", n_cuts, -0.5, n_cuts - 0.5)
    hist_ratio_sps_jpsiincl_dmesincl = ROOT.TH1D("hist_ratio_sps_jpsiincl_dmesincl", ";;#it{R}_{DJ/#Psi} (#mub)", n_cuts, -0.5, n_cuts - 0.5)
    hist_ratio_dps_jpsiincl_dmesincl = ROOT.TH1D("hist_ratio_dps_jpsiincl_dmesincl", ";;#it{R}_{DJ/#Psi} (#mub)", n_cuts, -0.5, n_cuts - 0.5)
    hist_ratio_tot_jpsiincl_dmesprompt = ROOT.TH1D("hist_ratio_tot_jpsiincl_dmesprompt", ";;#it{R}_{DJ/#Psi} (#mub)", n_cuts, -0.5, n_cuts-0.5)
    hist_ratio_sps_jpsiincl_dmesprompt = ROOT.TH1D("hist_ratio_sps_jpsiincl_dmesprompt", ";;#it{R}_{DJ/#Psi} (#mub)", n_cuts, -0.5, n_cuts-0.5)
    hist_ratio_dps_jpsiincl_dmesprompt = ROOT.TH1D("hist_ratio_dps_jpsiincl_dmesprompt", ";;#it{R}_{DJ/#Psi} (#mub)", n_cuts, -0.5, n_cuts-0.5)
    hist_ratio_tot_jpsiprompt_dmesprompt = ROOT.TH1D("hist_ratio_tot_jpsiprompt_dmesprompt", ";;#it{R}_{DJ/#Psi} (#mub)", n_cuts, -0.5, n_cuts - 0.5)
    hist_ratio_sps_jpsiprompt_dmesprompt = ROOT.TH1D("hist_ratio_sps_jpsiprompt_dmesprompt", ";;#it{R}_{DJ/#Psi} (#mub)", n_cuts, -0.5, n_cuts - 0.5)
    hist_ratio_dps_jpsiprompt_dmesprompt = ROOT.TH1D("hist_ratio_dps_jpsiprompt_dmesprompt", ";;#it{R}_{DJ/#Psi} (#mub)", n_cuts, -0.5, n_cuts-0.5)
    set_style(hist_ratio_sps_jpsiincl_dmesincl, COLORS[0])
    set_style(hist_ratio_dps_jpsiincl_dmesincl, COLORS[1])
    set_style(hist_ratio_tot_jpsiincl_dmesincl, COLORS[2])
    set_style(hist_ratio_sps_jpsiincl_dmesprompt, COLORS[0])
    set_style(hist_ratio_dps_jpsiincl_dmesprompt, COLORS[1])
    set_style(hist_ratio_tot_jpsiincl_dmesprompt, COLORS[2])
    set_style(hist_ratio_sps_jpsiprompt_dmesprompt, COLORS[0])
    set_style(hist_ratio_dps_jpsiprompt_dmesprompt, COLORS[1])
    set_style(hist_ratio_tot_jpsiprompt_dmesprompt, COLORS[2])
    hist_ratio_sps_jpsiincl_dmesincl.GetXaxis().SetLabelSize(0.05)
    hist_ratio_dps_jpsiincl_dmesincl.GetXaxis().SetLabelSize(0.05)
    hist_ratio_tot_jpsiincl_dmesincl.GetXaxis().SetLabelSize(0.05)
    hist_ratio_sps_jpsiincl_dmesprompt.GetXaxis().SetLabelSize(0.05)
    hist_ratio_dps_jpsiincl_dmesprompt.GetXaxis().SetLabelSize(0.05)
    hist_ratio_tot_jpsiincl_dmesprompt.GetXaxis().SetLabelSize(0.05)
    hist_ratio_sps_jpsiprompt_dmesprompt.GetXaxis().SetLabelSize(0.05)
    hist_ratio_dps_jpsiprompt_dmesprompt.GetXaxis().SetLabelSize(0.05)
    hist_ratio_tot_jpsiprompt_dmesprompt.GetXaxis().SetLabelSize(0.05)

    if plot_lhcb and "lhcb_cuts" not in kine_cuts:
        print("WARNING: disabling comparison with LHCb, since no corresponding kine cuts enabled")
        plot_lhcb = False

    for ivar, variant in enumerate(cuts):

        hist_ratio_tot_jpsiincl_dmesincl.GetXaxis().SetBinLabel(ivar+1, variant)
        hist_ratio_sps_jpsiincl_dmesincl.GetXaxis().SetBinLabel(ivar+1, variant)
        hist_ratio_dps_jpsiincl_dmesincl.GetXaxis().SetBinLabel(ivar+1, variant)
        hist_ratio_tot_jpsiincl_dmesprompt.GetXaxis().SetBinLabel(ivar+1, variant)
        hist_ratio_sps_jpsiincl_dmesprompt.GetXaxis().SetBinLabel(ivar+1, variant)
        hist_ratio_dps_jpsiincl_dmesprompt.GetXaxis().SetBinLabel(ivar+1, variant)
        hist_ratio_tot_jpsiprompt_dmesprompt.GetXaxis().SetBinLabel(ivar+1, variant)
        hist_ratio_sps_jpsiprompt_dmesprompt.GetXaxis().SetBinLabel(ivar+1, variant)
        hist_ratio_dps_jpsiprompt_dmesprompt.GetXaxis().SetBinLabel(ivar+1, variant)

        canv[variant] = ROOT.TCanvas(f"canv_{variant}", "", 1000, 500)
        canv[variant].Divide(2, 1)

        header = "D^{0}#minusJ/#Psi (no kine cuts)"
        if variant == "fwd_cuts":
            header = "D^{0}(mid)-J/#Psi(fwd)"
        elif variant == "mid_cuts":
            header = "D^{0}(mid)-J/#Psi(mid)"
        elif variant == "lhcb_cuts":
            header = "D^{0}-J/#Psi LHCb"
        leg[variant] = ROOT.TLegend(0.2, 0.7, 0.5, 0.9)
        leg[variant].SetFillStyle(0)
        leg[variant].SetBorderSize(0)
        leg[variant].SetTextSize(0.045)
        leg[variant].SetHeader(header)
        leg[variant].AddEntry(hist_ratio_tot_jpsiincl_dmesincl, "Total", "lp")
        leg[variant].AddEntry(hist_ratio_sps_jpsiincl_dmesincl, "SPS", "lp")
        leg[variant].AddEntry(hist_ratio_dps_jpsiincl_dmesincl, "DPS", "lp")

        cuts_string = ""
        xsec_jpsi, xsec_dmes, xsec_dmes_prompt, xsec_jpsi_prompt = 0., 0., 0., 0.
        ptjpsi_min, ptjpsi_max, yjpsi_min, yjpsi_max = -1, -1, -1, -1
        ptdmes_min, ptdmes_max, ydmes_min, ydmes_max = -1, -1, -1, -1
        if variant == "no_kine_cuts":
            yjpsi_min = 1
            yjpsi_max = hist_jpsi.GetYaxis().GetNbins() + 1
            ptjpsi_min = 1
            ptjpsi_max = hist_jpsi.GetXaxis().GetNbins() + 1
            ydmes_min = 1
            ydmes_max = hist_dmes.GetYaxis().GetNbins() + 1
            ptdmes_min = 1
            ptdmes_max = hist_dmes.GetXaxis().GetNbins() + 1
        elif variant == "fwd_cuts":
            cuts_string = "and 0.5 < pTD < 24 and abs(yD) < 0.6 and pTJPsi < 20 and 2.5 < yJPsi < 4.0"
            yjpsi_min = hist_jpsi.GetYaxis().FindBin(2.5001)
            yjpsi_max = hist_jpsi.GetYaxis().FindBin(3.9999)
            ptjpsi_min = 1
            ptjpsi_max = hist_jpsi.GetXaxis().FindBin(19.9999)
            ydmes_min = hist_dmes.GetYaxis().FindBin(-0.5999)
            ydmes_max = hist_dmes.GetYaxis().FindBin(0.5999)
            ptdmes_min = hist_dmes.GetXaxis().FindBin(0.5001)
            ptdmes_max = hist_dmes.GetXaxis().FindBin(23.9999)
        elif variant == "mid_cuts":
            cuts_string = "and 0.5 < pTD < 24 and abs(yD) < 0.6 and pTJPsi < 20 and abs(yJPsi) < 0.9"
            yjpsi_min = hist_jpsi.GetYaxis().FindBin(-0.8999)
            yjpsi_max = hist_jpsi.GetYaxis().FindBin(0.8999)
            ptjpsi_min = 1
            ptjpsi_max = hist_jpsi.GetXaxis().FindBin(19.9999)
            ydmes_min = hist_dmes.GetYaxis().FindBin(-0.5999)
            ydmes_max = hist_dmes.GetYaxis().FindBin(0.5999)
            ptdmes_min = hist_dmes.GetXaxis().FindBin(0.5001)
            ptdmes_max = hist_dmes.GetXaxis().FindBin(23.9999)
        elif variant == "lhcb_cuts":
            cuts_string = "and 3 < pTD < 12 and 2 < yD < 4 and pTJPsi < 12 and 2 < yJPsi < 4"
            yjpsi_min = hist_jpsi.GetYaxis().FindBin(2.0001)
            yjpsi_max = hist_jpsi.GetYaxis().FindBin(4.0001)
            ptjpsi_min = 1
            ptjpsi_max = hist_jpsi.GetXaxis().FindBin(11.9999)
            ydmes_min = hist_dmes.GetYaxis().FindBin(2.0001)
            ydmes_max = hist_dmes.GetYaxis().FindBin(4.0001)
            ptdmes_min = hist_dmes.GetXaxis().FindBin(3.0001)
            ptdmes_max = hist_dmes.GetXaxis().FindBin(11.9999)

        xsec_jpsi = hist_jpsi.Integral(ptjpsi_min, ptjpsi_max, yjpsi_min, yjpsi_max, 1, 2) * xsec / nevents
        xsec_jpsi_prompt = hist_jpsi.Integral(ptjpsi_min, ptjpsi_max, yjpsi_min, yjpsi_max, 1, 1) * xsec / nevents
        xsec_dmes = hist_dmes.Integral(ptdmes_min, ptdmes_max, ydmes_min, ydmes_max, 1, 2) * xsec / nevents
        xsec_dmes_prompt = hist_dmes.Integral(ptdmes_min, ptdmes_max, ydmes_min, ydmes_max, 1, 1) * xsec / nevents

        df["deltaY"] = df["yJPsi"] - df["yD"]
        df["deltaPhi"] = df.apply(lambda row: abs(row.phiJPsi - row.phiD) / np.pi if abs(row.phiJPsi - row.phiD) < np.pi else (2*np.pi - abs(row.phiJPsi - row.phiD)) / np.pi, axis=1)

        df_sps = df.query(f"commonAnchestor > 0.5 and abs(pdgD) == 421{cuts_string}")
        hist_deltay_sps_jpsiincl_dmesincl = ROOT.TH1D("hist_deltay_sps_jpsiincl_dmesincl",
                                                      ";#Delta#it{y}; d#sigma/d#Delta#it{y} (#mub)",
                                                      80, -10., 10.)
        hist_deltay_sps_jpsiincl_dmesprompt = ROOT.TH1D("hist_deltay_sps_jpsiincl_dmesprompt",
                                                        ";#Delta#it{y}; d#sigma/d#Delta#it{y} (#mub)",
                                                        80, -10., 10.)
        hist_deltay_sps_jpsiprompt_dmesprompt = ROOT.TH1D("hist_deltay_sps_jpsiprompt_dmesprompt",
                                                          ";#Delta#it{y}; d#sigma/d#Delta#it{y} (#mub)",
                                                          80, -10., 10.)
        hist_deltaphi_sps_jpsiincl_dmesincl = ROOT.TH1D("hist_deltaphi_sps_jpsiincl_dmesincl",
                                                        ";|#Delta#varphi|/#pi; d#sigma/(d|#Delta#varphi|/#pi) (#mub)",
                                                        10, 0., 1)
        hist_deltaphi_sps_jpsiincl_dmesprompt = ROOT.TH1D("hist_deltaphi_sps_jpsiincl_dmesprompt",
                                                          ";|#Delta#varphi|/#pi; d#sigma/|d#Delta#varphi| (#mub)",
                                                          10, 0., 1)
        hist_deltaphi_sps_jpsiprompt_dmesprompt = ROOT.TH1D("hist_deltaphi_sps_jpsiprompt_dmesprompt",
                                                            ";|#Delta#varphi|/#pi; d#sigma/|d#Delta#varphi| (#mub)",
                                                            10, 0., 1)

        for deltay, deltaphi, dfromb, jpsifromb in zip(df_sps["deltaY"].to_numpy(),
                                                       df_sps["deltaPhi"].to_numpy(),
                                                       df_sps["DfromB"].to_numpy(),
                                                       df_sps["JPsifromB"].to_numpy()):
            hist_deltay_sps_jpsiincl_dmesincl.Fill(deltay)
            hist_deltaphi_sps_jpsiincl_dmesincl.Fill(deltaphi)
            if not dfromb:
                hist_deltay_sps_jpsiincl_dmesprompt.Fill(deltay)
                hist_deltaphi_sps_jpsiincl_dmesprompt.Fill(deltaphi)
                if not jpsifromb:
                    hist_deltay_sps_jpsiprompt_dmesprompt.Fill(deltay)
                    hist_deltaphi_sps_jpsiprompt_dmesprompt.Fill(deltaphi)
        
        xsec_sps_jpsidmesincl = hist_deltay_sps_jpsiincl_dmesincl.Integral() * xsec / nevents / 2 # / 2 is for charge average
        unc_sps_incl = np.sqrt(hist_deltay_sps_jpsiincl_dmesincl.Integral()) * xsec / nevents / 2
        xsec_sps_jpsiincl_dmesprompt = hist_deltay_sps_jpsiincl_dmesprompt.Integral() * xsec / nevents / 2
        unc_sps_dprompt = np.sqrt(hist_deltay_sps_jpsiincl_dmesprompt.Integral()) * xsec / nevents / 2
        xsec_sps_jpsidmesprompt = hist_deltay_sps_jpsiprompt_dmesprompt.Integral() * xsec / nevents / 2
        unc_sps_prompt = np.sqrt(hist_deltay_sps_jpsiprompt_dmesprompt.Integral()) * xsec / nevents / 2
        ratio_sps_incl = xsec_jpsi * xsec_dmes / xsec_sps_jpsidmesincl
        ratiounc_sps_incl = unc_sps_incl / xsec_sps_jpsidmesincl * ratio_sps_incl
        hist_ratio_sps_jpsiincl_dmesincl.SetBinContent(ivar+1, ratio_sps_incl)
        hist_ratio_sps_jpsiincl_dmesincl.SetBinError(ivar+1, ratiounc_sps_incl)
        ratio_sps_dprompt = xsec_jpsi * xsec_dmes_prompt / xsec_sps_jpsiincl_dmesprompt
        ratiounc_sps_dprompt = unc_sps_dprompt / xsec_sps_jpsiincl_dmesprompt * ratio_sps_dprompt
        hist_ratio_sps_jpsiincl_dmesprompt.SetBinContent(ivar+1, ratio_sps_dprompt)
        hist_ratio_sps_jpsiincl_dmesprompt.SetBinError(ivar+1, ratiounc_sps_dprompt)
        ratio_sps_prompt = xsec_jpsi_prompt * xsec_dmes_prompt / xsec_sps_jpsidmesprompt
        ratiounc_sps_prompt = unc_sps_prompt / xsec_sps_jpsidmesprompt * ratio_sps_dprompt
        hist_ratio_sps_jpsiprompt_dmesprompt.SetBinContent(ivar+1, ratio_sps_prompt)
        hist_ratio_sps_jpsiprompt_dmesprompt.SetBinError(ivar+1, ratiounc_sps_prompt)

        hist_deltay_sps_jpsiincl_dmesincl.Scale(xsec/nevents, "width")
        hist_deltay_sps_jpsiincl_dmesprompt.Scale(xsec/nevents, "width")
        hist_deltay_sps_jpsiprompt_dmesprompt.Scale(xsec/nevents, "width")
        hist_deltaphi_sps_jpsiincl_dmesincl.Scale(xsec/nevents, "width")
        hist_deltaphi_sps_jpsiincl_dmesprompt.Scale(xsec/nevents, "width")
        hist_deltaphi_sps_jpsiprompt_dmesprompt.Scale(xsec/nevents, "width")

        set_style(hist_deltay_sps_jpsiincl_dmesincl, COLORS[0])
        set_style(hist_deltay_sps_jpsiincl_dmesprompt, COLORS[0])
        set_style(hist_deltay_sps_jpsiprompt_dmesprompt, COLORS[0])
        set_style(hist_deltaphi_sps_jpsiincl_dmesincl, COLORS[0])
        set_style(hist_deltaphi_sps_jpsiincl_dmesprompt, COLORS[0])
        set_style(hist_deltaphi_sps_jpsiprompt_dmesprompt, COLORS[0])

        df_dps = df.query(f"commonAnchestor < 0.5 and abs(pdgD) == 421{cuts_string}")
        hist_deltay_dps_jpsiincl_dmesincl = ROOT.TH1D("hist_deltay_dps_jpsiincl_dmesincl",
                                                      ";#Delta#it{y}; d#sigma/d#Delta#it{y} (#mub)",
                                                      80, -10., 10.)
        hist_deltay_dps_jpsiincl_dmesprompt = ROOT.TH1D("hist_deltay_dps_jpsiincl_dmesprompt",
                                                        ";#Delta#it{y}; d#sigma/d#Delta#it{y} (#mub)",
                                                        80, -10., 10.)
        hist_deltay_dps_jpsiprompt_dmesprompt = ROOT.TH1D("hist_deltay_dps_jpsiprompt_dmesprompt",
                                                          ";#Delta#it{y}; d#sigma/d#Delta#it{y} (#mub)",
                                                          80, -10., 10.)
        hist_deltaphi_dps_jpsiincl_dmesincl = ROOT.TH1D("hist_deltaphi_dps_jpsiincl_dmesincl",
                                                        ";|#Delta#varphi|/#pi; d#sigma/(d|#Delta#varphi|/#pi) (#mub)",
                                                        10, 0., 1.)
        hist_deltaphi_dps_jpsiincl_dmesprompt = ROOT.TH1D("hist_deltaphi_dps_jpsiincl_dmesprompt",
                                                          ";|#Delta#varphi|/#pi; d#sigma/(d|#Delta#varphi|/#pi) (#mub)",
                                                          10, 0., 1.)
        hist_deltaphi_dps_jpsiprompt_dmesprompt = ROOT.TH1D("hist_deltaphi_dps_jpsiprompt_dmesprompt",
                                                            ";|#Delta#varphi|/#pi; d#sigma/|d#Delta#varphi| (#mub)",
                                                            10, 0., 1.)
        for deltay, deltaphi, dfromb, jpsifromb in zip(df_dps["deltaY"].to_numpy(),
                                                      df_dps["deltaPhi"].to_numpy(),
                                                      df_dps["DfromB"].to_numpy(),
                                                      df_dps["JPsifromB"].to_numpy()):
            hist_deltay_dps_jpsiincl_dmesincl.Fill(deltay)
            hist_deltaphi_dps_jpsiincl_dmesincl.Fill(deltaphi)
            if not dfromb:
                hist_deltay_dps_jpsiincl_dmesprompt.Fill(deltay)
                hist_deltaphi_dps_jpsiincl_dmesprompt.Fill(deltaphi)
                if not jpsifromb:
                    hist_deltay_dps_jpsiprompt_dmesprompt.Fill(deltay)
                    hist_deltaphi_dps_jpsiprompt_dmesprompt.Fill(deltaphi)

        xsec_dps_jpsidmesincl = hist_deltay_dps_jpsiincl_dmesincl.Integral() * xsec / nevents / 2
        unc_dps_incl = np.sqrt(hist_deltay_dps_jpsiincl_dmesincl.Integral()) * xsec / nevents / 2
        xsec_dps_jpsiincl_dmesprompt = hist_deltay_dps_jpsiincl_dmesprompt.Integral() * xsec / nevents / 2
        unc_dps_dprompt = np.sqrt(hist_deltay_dps_jpsiincl_dmesprompt.Integral()) * xsec / nevents / 2
        xsec_dps_jpsidmesprompt = hist_deltay_dps_jpsiprompt_dmesprompt.Integral() * xsec / nevents / 2
        unc_dps_prompt = np.sqrt(hist_deltay_dps_jpsiprompt_dmesprompt.Integral()) * xsec / nevents / 2
        ratio_dps_incl = xsec_jpsi * xsec_dmes / xsec_dps_jpsidmesincl
        ratio_dps_dprompt = xsec_jpsi * xsec_dmes_prompt / xsec_dps_jpsiincl_dmesprompt
        ratio_dps_prompt = xsec_jpsi_prompt * xsec_dmes_prompt / xsec_dps_jpsidmesprompt
        ratiounc_dps_incl = unc_dps_incl / xsec_dps_jpsidmesincl * ratio_dps_incl
        ratiounc_dps_dprompt = unc_dps_dprompt / xsec_dps_jpsiincl_dmesprompt * ratio_dps_dprompt
        ratiounc_dps_prompt = unc_dps_prompt / xsec_dps_jpsidmesprompt * ratio_dps_prompt
        unc_tot_incl = np.sqrt(hist_deltay_dps_jpsiincl_dmesincl.Integral() + xsec_sps_jpsidmesincl * nevents / xsec) * xsec / nevents
        unc_tot_dprompt = np.sqrt(hist_deltay_dps_jpsiincl_dmesprompt.Integral() + xsec_sps_jpsiincl_dmesprompt * nevents / xsec) * xsec / nevents
        unc_tot_prompt = np.sqrt(hist_deltay_dps_jpsiprompt_dmesprompt.Integral() + xsec_sps_jpsidmesprompt * nevents / xsec) * xsec / nevents
        ratio_tot_incl = xsec_jpsi * xsec_dmes / (xsec_sps_jpsidmesincl + xsec_dps_jpsidmesincl)
        ratio_tot_dprompt = xsec_jpsi * xsec_dmes_prompt / (xsec_sps_jpsiincl_dmesprompt + xsec_dps_jpsiincl_dmesprompt)
        ratio_tot_prompt = xsec_jpsi_prompt * xsec_dmes_prompt / (xsec_sps_jpsidmesprompt + xsec_dps_jpsidmesprompt)
        ratiounc_tot_incl = unc_tot_incl / (xsec_sps_jpsidmesincl + xsec_dps_jpsidmesincl) * ratio_tot_incl
        ratiounc_tot_dprompt = unc_tot_dprompt / (xsec_sps_jpsiincl_dmesprompt + xsec_dps_jpsiincl_dmesprompt) * ratio_tot_dprompt
        ratiounc_tot_prompt = unc_tot_prompt / (xsec_sps_jpsidmesprompt + xsec_dps_jpsidmesprompt) * ratio_tot_prompt
        hist_ratio_dps_jpsiincl_dmesincl.SetBinContent(ivar+1, ratio_dps_incl)
        hist_ratio_dps_jpsiincl_dmesprompt.SetBinContent(ivar+1, ratio_dps_dprompt)
        hist_ratio_dps_jpsiprompt_dmesprompt.SetBinContent(ivar+1, ratio_dps_prompt)
        hist_ratio_tot_jpsiincl_dmesincl.SetBinContent(ivar+1, ratio_tot_incl)
        hist_ratio_tot_jpsiincl_dmesprompt.SetBinContent(ivar+1, ratio_tot_dprompt)
        hist_ratio_tot_jpsiprompt_dmesprompt.SetBinContent(ivar+1, ratio_tot_prompt)
        hist_ratio_dps_jpsiincl_dmesincl.SetBinError(ivar+1, ratiounc_dps_incl)
        hist_ratio_dps_jpsiincl_dmesprompt.SetBinError(ivar+1, ratiounc_dps_dprompt)
        hist_ratio_dps_jpsiprompt_dmesprompt.SetBinError(ivar+1, ratiounc_dps_prompt)
        hist_ratio_tot_jpsiincl_dmesincl.SetBinError(ivar+1, ratiounc_tot_incl)
        hist_ratio_tot_jpsiincl_dmesprompt.SetBinError(ivar+1, ratiounc_tot_dprompt)
        hist_ratio_tot_jpsiprompt_dmesprompt.SetBinError(ivar+1, ratiounc_tot_prompt)

        hist_deltay_dps_jpsiincl_dmesincl.Scale(xsec/nevents, "width")
        hist_deltay_dps_jpsiincl_dmesprompt.Scale(xsec/nevents, "width")
        hist_deltay_dps_jpsiprompt_dmesprompt.Scale(xsec/nevents, "width")
        hist_deltaphi_dps_jpsiincl_dmesincl.Scale(xsec/nevents, "width")
        hist_deltaphi_dps_jpsiincl_dmesprompt.Scale(xsec/nevents, "width")
        hist_deltaphi_dps_jpsiprompt_dmesprompt.Scale(xsec/nevents, "width")

        set_style(hist_deltay_dps_jpsiincl_dmesincl, COLORS[1])
        set_style(hist_deltay_dps_jpsiincl_dmesprompt, COLORS[1])
        set_style(hist_deltay_dps_jpsiprompt_dmesprompt, COLORS[1])
        set_style(hist_deltaphi_dps_jpsiincl_dmesincl, COLORS[1])
        set_style(hist_deltaphi_dps_jpsiincl_dmesprompt, COLORS[1])
        set_style(hist_deltaphi_dps_jpsiprompt_dmesprompt, COLORS[1])

        hist_deltay_tot_jpsiincl_dmesincl = hist_deltay_dps_jpsiincl_dmesincl.Clone(
            "hist_deltay_tot_jpsiincl_dmesincl")
        hist_deltay_tot_jpsiincl_dmesprompt = hist_deltay_dps_jpsiincl_dmesprompt.Clone(
            "hist_deltay_tot_jpsiincl_dmesprompt")
        hist_deltay_tot_jpsiprompt_dmesprompt = hist_deltay_dps_jpsiprompt_dmesprompt.Clone(
            "hist_deltay_tot_jpsiprompt_dmesprompt")
        hist_deltaphi_tot_jpsiincl_dmesincl = hist_deltaphi_dps_jpsiincl_dmesincl.Clone(
            "hist_deltaphi_tot_jpsiincl_dmesincl")
        hist_deltaphi_tot_jpsiincl_dmesprompt = hist_deltaphi_dps_jpsiincl_dmesprompt.Clone(
            "hist_deltaphi_tot_jpsiincl_dmesprompt")
        hist_deltaphi_tot_jpsiprompt_dmesprompt = hist_deltaphi_dps_jpsiprompt_dmesprompt.Clone(
            "hist_deltaphi_tot_jpsiprompt_dmesprompt")
        hist_deltay_tot_jpsiincl_dmesincl.Add(hist_deltay_sps_jpsiincl_dmesincl)
        hist_deltay_tot_jpsiincl_dmesprompt.Add(hist_deltay_sps_jpsiincl_dmesprompt)
        hist_deltay_tot_jpsiprompt_dmesprompt.Add(hist_deltay_sps_jpsiprompt_dmesprompt)
        hist_deltaphi_tot_jpsiincl_dmesincl.Add(hist_deltaphi_sps_jpsiincl_dmesincl)
        hist_deltaphi_tot_jpsiincl_dmesprompt.Add(hist_deltaphi_sps_jpsiincl_dmesprompt)
        hist_deltaphi_tot_jpsiprompt_dmesprompt.Add(hist_deltaphi_sps_jpsiprompt_dmesprompt)

        set_style(hist_deltay_tot_jpsiincl_dmesincl, COLORS[2])
        set_style(hist_deltay_tot_jpsiincl_dmesprompt, COLORS[2])
        set_style(hist_deltay_tot_jpsiprompt_dmesprompt, COLORS[2])
        set_style(hist_deltaphi_tot_jpsiincl_dmesincl, COLORS[2])
        set_style(hist_deltaphi_tot_jpsiincl_dmesprompt, COLORS[2])
        set_style(hist_deltaphi_tot_jpsiprompt_dmesprompt, COLORS[2])

        canv[variant].cd(1)
        hist_deltay_tot_jpsiincl_dmesprompt.GetYaxis().SetDecimals()
        hist_deltay_tot_jpsiincl_dmesprompt.GetYaxis().SetRangeUser(
            0., hist_deltay_tot_jpsiincl_dmesprompt.GetMaximum()*1.5)
        hist_deltay_tot_jpsiincl_dmesprompt.DrawCopy()
        hist_deltay_sps_jpsiincl_dmesprompt.DrawCopy("same")
        hist_deltay_dps_jpsiincl_dmesprompt.DrawCopy("same")
        leg[variant].Draw()
        canv[variant].cd(2)
        hist_deltaphi_tot_jpsiincl_dmesprompt.GetYaxis().SetDecimals()
        hist_deltaphi_tot_jpsiincl_dmesprompt.GetYaxis().SetRangeUser(
            0., hist_deltaphi_tot_jpsiincl_dmesprompt.GetMaximum()*1.2)
        hist_deltaphi_tot_jpsiincl_dmesprompt.DrawCopy()
        hist_deltaphi_sps_jpsiincl_dmesprompt.DrawCopy("same")
        hist_deltaphi_dps_jpsiincl_dmesprompt.DrawCopy("same")

        canv[variant].SaveAs(outfile_name.replace(".root", f"_{variant}.pdf"))

        if plot_lhcb and variant == "lhcb_cuts":
            lhcb_data = ROOT.TFile.Open("HEPData-ins1113596-v1-root.root")
            # deltaphi distribution
            graph_deltaphi_lhcb = lhcb_data.Get("Table 21/Graph1D_y1")
            graph_deltay_lhcb = lhcb_data.Get("Table 24/Graph1D_y1")
            graph_deltaphi_lhcb.SetNameTitle("graph_deltaphi_lhcb",
                                             ";|#Delta#varphi|/#pi;d#sigma/(d|#Delta#varphi|/#pi) (a.u.)")
            graph_deltay_lhcb.SetNameTitle("graph_deltay_lhcb", ";#Delta#it{y};d#sigma/(d#Delta#it{y}) (a.u.)")
            ratio_lhcb = 14.9e3
            unc_tot_low_ratio_lhcb = np.sqrt(0.4**2 + 1.1**2 + 3.1**2)*1.e3
            unc_tot_high_ratio_lhcb = np.sqrt(0.4**2 + 1.1**2 + 2.3**2)*1.e3
            graph_ratio_lhcb = ROOT.TGraphAsymmErrors(1)
            graph_ratio_lhcb.SetNameTitle("graph_ratio_lhcb", ";;#it{R}_{DJ/#Psi}")
            graph_ratio_lhcb.SetPoint(0, 1., ratio_lhcb)
            graph_ratio_lhcb.SetPointError(0, 0., 0., unc_tot_low_ratio_lhcb, unc_tot_high_ratio_lhcb)
            set_style(graph_deltaphi_lhcb, ROOT.kBlack, ROOT.kOpenSquare, 1.5)
            set_style(graph_deltay_lhcb, ROOT.kBlack, ROOT.kOpenSquare, 1.5)
            set_style(graph_ratio_lhcb, ROOT.kBlack, ROOT.kOpenSquare, 1.5)
            canv_lhcb = ROOT.TCanvas("canv_lhcb", "", 1500, 500)
            canv_lhcb.Divide(3, 1)
            frame = canv_lhcb.cd(1).DrawFrame(0., 0., 1., 0.25, ";|#Delta#varphi|/#pi;d#sigma/(d|#Delta#varphi|/#pi) (a.u.)")
            frame.GetYaxis().SetDecimals()
            hist_deltaphi_tot_pythia4lhcb = hist_deltaphi_tot_jpsiincl_dmesprompt.Clone("hist_deltaphi_tot_pythia4lhcb")
            hist_deltaphi_sps_pythia4lhcb = hist_deltaphi_sps_jpsiincl_dmesprompt.Clone("hist_deltaphi_sps_pythia4lhcb")
            hist_deltaphi_dps_pythia4lhcb = hist_deltaphi_dps_jpsiincl_dmesprompt.Clone("hist_deltaphi_dps_pythia4lhcb")
            norm = 1. / hist_deltaphi_tot_pythia4lhcb.Integral()
            hist_deltaphi_tot_pythia4lhcb.Scale(norm)
            hist_deltaphi_sps_pythia4lhcb.Scale(norm)
            hist_deltaphi_dps_pythia4lhcb.Scale(norm)
            hist_deltaphi_tot_pythia4lhcb.DrawCopy("same")
            hist_deltaphi_sps_pythia4lhcb.DrawCopy("same")
            hist_deltaphi_dps_pythia4lhcb.DrawCopy("same")
            graph_deltaphi_lhcb.Draw("pz")
            leg_lhcb_pythia = ROOT.TLegend(0.2, 0.6, 0.5, 0.85)
            leg_lhcb_pythia.SetBorderSize(0)
            leg_lhcb_pythia.SetFillStyle(0)
            leg_lhcb_pythia.SetTextSize(0.045)
            leg_lhcb_pythia.SetHeader("PYTHIA 8, D^{0}-J/#Psi")
            leg_lhcb_pythia.AddEntry(hist_deltaphi_tot_pythia4lhcb, "Total", "lp")
            leg_lhcb_pythia.AddEntry(hist_deltaphi_sps_pythia4lhcb, "SPS", "lp")
            leg_lhcb_pythia.AddEntry(hist_deltaphi_dps_pythia4lhcb, "DPS", "lp")
            leg_lhcb_pythia.Draw()
            frame = canv_lhcb.cd(2).DrawFrame(-2., 0., 2., 0.25, ";#Delta#it{y};d#sigma/(d#Delta#it{y}) (a.u.)")
            frame.GetYaxis().SetDecimals()
            hist_deltay_tot_pythia4lhcb = hist_deltay_tot_jpsiincl_dmesprompt.Clone("hist_deltay_tot_pythia4lhcb")
            hist_deltay_sps_pythia4lhcb = hist_deltay_sps_jpsiincl_dmesprompt.Clone("hist_deltay_sps_pythia4lhcb")
            hist_deltay_dps_pythia4lhcb = hist_deltay_dps_jpsiincl_dmesprompt.Clone("hist_deltay_dps_pythia4lhcb")
            norm = 1. / hist_deltay_tot_pythia4lhcb.Integral()
            hist_deltay_tot_pythia4lhcb.Scale(norm)
            hist_deltay_sps_pythia4lhcb.Scale(norm)
            hist_deltay_dps_pythia4lhcb.Scale(norm)
            hist_deltay_tot_pythia4lhcb.DrawCopy("same")
            hist_deltay_sps_pythia4lhcb.DrawCopy("same")
            hist_deltay_dps_pythia4lhcb.DrawCopy("same")
            graph_deltay_lhcb.Draw("pz")
            leg_lhcb_data = ROOT.TLegend(0.2, 0.75, 0.5, 0.85)
            leg_lhcb_data.SetBorderSize(0)
            leg_lhcb_data.SetFillStyle(0)
            leg_lhcb_data.SetTextSize(0.045)
            leg_lhcb_data.AddEntry(graph_deltay_lhcb, "LHCb D^{0}-J/#Psi", "lp")
            leg_lhcb_data.AddEntry("", "LHCb JHEP 06 (2012) 141", "")
            leg_lhcb_data.Draw()
            canv_lhcb.cd(3).SetLogy()
            hist_ratio_dps_jpsiprompt_dmesprompt.GetYaxis().SetRangeUser(1.e3, 4.e5)
            hist_ratio_dps_jpsiprompt_dmesprompt.DrawCopy()
            hist_ratio_sps_jpsiprompt_dmesprompt.DrawCopy("same")
            hist_ratio_tot_jpsiprompt_dmesprompt.DrawCopy("same")
            graph_ratio_lhcb.Draw("pz")

            canv_lhcb.SaveAs(outfile_name.replace(".root", "_comparison_LHCb.pdf"))

            outfile.cd()
            outfile.mkdir("lhcb_comparison")
            canv_lhcb.Write()
            hist_deltaphi_tot_pythia4lhcb.Write()
            hist_deltaphi_sps_pythia4lhcb.Write()
            hist_deltaphi_dps_pythia4lhcb.Write()
            hist_deltay_tot_pythia4lhcb.Write()
            hist_deltay_sps_pythia4lhcb.Write()
            hist_deltay_dps_pythia4lhcb.Write()
            hist_ratio_tot_jpsiincl_dmesprompt.Write()
            graph_deltaphi_lhcb
            graph_deltay_lhcb.Write()
            graph_ratio_lhcb.Write()

        outfile.cd()
        outfile.mkdir(variant)
        outfile.cd(variant)
        canv[variant].Write()
        hist_deltay_sps_jpsiincl_dmesincl.Write()
        hist_deltay_sps_jpsiincl_dmesprompt.Write()
        hist_deltay_sps_jpsiprompt_dmesprompt.Write()
        hist_deltaphi_sps_jpsiincl_dmesincl.Write()
        hist_deltaphi_sps_jpsiincl_dmesprompt.Write()
        hist_deltaphi_sps_jpsiprompt_dmesprompt.Write()
        hist_deltay_dps_jpsiincl_dmesincl.Write()
        hist_deltay_dps_jpsiincl_dmesprompt.Write()
        hist_deltay_dps_jpsiprompt_dmesprompt.Write()
        hist_deltaphi_dps_jpsiincl_dmesincl.Write()
        hist_deltaphi_dps_jpsiincl_dmesprompt.Write()
        hist_deltaphi_dps_jpsiprompt_dmesprompt.Write()
        hist_deltay_tot_jpsiincl_dmesincl.Write()
        hist_deltay_tot_jpsiincl_dmesprompt.Write()
        hist_deltay_tot_jpsiprompt_dmesprompt.Write()
        hist_deltaphi_tot_jpsiincl_dmesincl.Write()
        hist_deltaphi_tot_jpsiincl_dmesprompt.Write()
        hist_deltaphi_tot_jpsiprompt_dmesprompt.Write()

    canv_ratio = ROOT.TCanvas("canv_ratio", "", 1000, 500)
    canv_ratio.SetLogy()
    hist_ratio_dps_jpsiincl_dmesprompt.GetYaxis().SetRangeUser(1.e2, 1.e6)
    hist_ratio_dps_jpsiincl_dmesprompt.DrawCopy()
    hist_ratio_sps_jpsiincl_dmesprompt.DrawCopy("same")
    hist_ratio_tot_jpsiincl_dmesprompt.DrawCopy("same")
    canv_ratio.SaveAs(outfile_name.replace(".root", "_ratio.pdf"))

    outfile.cd()
    canv_ratio.Write()
    hist_ratio_dps_jpsiincl_dmesincl.Write()
    hist_ratio_sps_jpsiincl_dmesincl.Write()
    hist_ratio_tot_jpsiincl_dmesincl.Write()
    hist_ratio_dps_jpsiincl_dmesprompt.Write()
    hist_ratio_sps_jpsiincl_dmesprompt.Write()
    hist_ratio_tot_jpsiincl_dmesprompt.Write()
    outfile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("-i", "--infile", metavar="text", required=True,
                        help="ROOT input file")
    parser.add_argument("-o", "--outfile", metavar="text", required=True,
                        help="ROOT output file name")
    parser.add_argument("-k", "--kine_cuts", nargs="+", type=str, help="list of kine cuts to be applied",
                        default=["fwd_cuts", "mid_cuts"], required=False) # options: ["fwd_cuts", "mid_cuts", "lhcb_cuts"]
    parser.add_argument("-p", "--plot_lhcb", action="store_true", help="option to plot predictions vs LHCb data", 
                        default=False, required=False)
    args = parser.parse_args()

    produce_predictions(args.infile, args.outfile, args.kine_cuts, args.plot_lhcb)
