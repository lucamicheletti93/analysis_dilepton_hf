"""
"""

import argparse
import uproot
import numpy as np
import ROOT

COLORS = [ROOT.kRed+1, ROOT.kAzure+4, ROOT.kBlack]

def set_style(hist, color):
    """
    """

    hist.SetLineWidth(2)
    hist.SetMarkerStyle(ROOT.kFullCircle)
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)


def produce_predictions(infile_name, outfile_name):
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
    infile_root.Close()

    df = uproot.open(infile_name)["tuplePairs"].arrays(library="pd")

    outfile = ROOT.TFile(outfile_name, "recreate")
    canv = {}
    hist_ratio_tot_jpsiincl_dmesincl = ROOT.TH1D("hist_ratio_tot_jpsiincl_dmesincl", ";;#it{R}_{DJ/#Psi} (#mub)", 3, -0.5, 2.5)
    hist_ratio_sps_jpsiincl_dmesincl = ROOT.TH1D("hist_ratio_sps_jpsiincl_dmesincl", ";;#it{R}_{DJ/#Psi} (#mub)", 3, -0.5, 2.5)
    hist_ratio_dps_jpsiincl_dmesincl = ROOT.TH1D("hist_ratio_dps_jpsiincl_dmesincl", ";;#it{R}_{DJ/#Psi} (#mub)", 3, -0.5, 2.5)
    hist_ratio_tot_jpsiincl_dmesprompt = ROOT.TH1D("hist_ratio_tot_jpsiincl_dmesprompt", ";;#it{R}_{DJ/#Psi} (#mub)", 3, -0.5, 2.5)
    hist_ratio_sps_jpsiincl_dmesprompt = ROOT.TH1D("hist_ratio_sps_jpsiincl_dmesprompt", ";;#it{R}_{DJ/#Psi} (#mub)", 3, -0.5, 2.5)
    hist_ratio_dps_jpsiincl_dmesprompt = ROOT.TH1D("hist_ratio_dps_jpsiincl_dmesprompt", ";;#it{R}_{DJ/#Psi} (#mub)", 3, -0.5, 2.5)
    set_style(hist_ratio_sps_jpsiincl_dmesincl, COLORS[0])
    set_style(hist_ratio_dps_jpsiincl_dmesincl, COLORS[1])
    set_style(hist_ratio_tot_jpsiincl_dmesincl, COLORS[2])
    set_style(hist_ratio_sps_jpsiincl_dmesprompt, COLORS[0])
    set_style(hist_ratio_dps_jpsiincl_dmesprompt, COLORS[1])
    set_style(hist_ratio_tot_jpsiincl_dmesprompt, COLORS[2])
    hist_ratio_sps_jpsiincl_dmesincl.GetXaxis().SetLabelSize(0.05)
    hist_ratio_dps_jpsiincl_dmesprompt.GetXaxis().SetLabelSize(0.05)

    for ivar, variant in enumerate(["no_kine_cuts", "fwd_cuts", "mid_cuts"]):

        hist_ratio_tot_jpsiincl_dmesincl.GetXaxis().SetBinLabel(ivar+1, variant)
        hist_ratio_sps_jpsiincl_dmesincl.GetXaxis().SetBinLabel(ivar+1, variant)
        hist_ratio_dps_jpsiincl_dmesincl.GetXaxis().SetBinLabel(ivar+1, variant)
        hist_ratio_tot_jpsiincl_dmesprompt.GetXaxis().SetBinLabel(ivar+1, variant)
        hist_ratio_sps_jpsiincl_dmesprompt.GetXaxis().SetBinLabel(ivar+1, variant)
        hist_ratio_dps_jpsiincl_dmesprompt.GetXaxis().SetBinLabel(ivar+1, variant)

        canv[variant] = ROOT.TCanvas(f"canv_{variant}", "", 1000, 500)
        canv[variant].Divide(2, 1)

        kine_cuts = ""
        if variant == "no_kine_cuts":
            xsec_jpsi = hist_jpsi.Integral() * xsec / nevents
            xsec_dmes = hist_dmes.Integral() * xsec / nevents
            xsec_dmes_prompt = hist_dmes.Integral(1, 1000000, 1, 10000000, 1, 1) * xsec / nevents
        elif variant == "fwd_cuts":
            kine_cuts = "and 0.5 < pTD < 24 and abs(yD) < 0.6 and pTJPsi < 20 and 2.5 < yJPsi < 4.0"
            yjpsi_min = hist_jpsi.GetYaxis().FindBin(2.5001)
            yjpsi_max = hist_jpsi.GetYaxis().FindBin(3.9999)
            ptjpsi_min = 1
            ptjpsi_max = hist_jpsi.GetXaxis().FindBin(19.9999)
            ydmes_min = hist_dmes.GetYaxis().FindBin(-0.5999)
            ydmes_max = hist_dmes.GetYaxis().FindBin(0.5999)
            ptdmes_min = hist_dmes.GetXaxis().FindBin(0.5001)
            ptdmes_max = hist_dmes.GetXaxis().FindBin(23.9999)
            xsec_jpsi = hist_jpsi.Integral(yjpsi_min, yjpsi_max, ptjpsi_min, ptjpsi_max, 1, 2) * xsec / nevents
            xsec_dmes = hist_dmes.Integral(ptdmes_min, ptdmes_max, ydmes_min, ydmes_max, 1, 2) * xsec / nevents
            xsec_dmes_prompt = hist_dmes.Integral(ptdmes_min, ptdmes_max, ydmes_min, ydmes_max, 1, 1) * xsec / nevents
        elif variant == "mid_cuts":
            kine_cuts = "and 0.5 < pTD < 24 and abs(yD) < 0.6 and pTJPsi < 20 and abs(yJPsi) < 0.9"
            yjpsi_min = hist_jpsi.GetYaxis().FindBin(-0.8999)
            yjpsi_max = hist_jpsi.GetYaxis().FindBin(0.8999)
            ptjpsi_min = 1
            ptjpsi_max = hist_jpsi.GetXaxis().FindBin(19.9999)
            ydmes_min = hist_dmes.GetYaxis().FindBin(-0.5999)
            ydmes_max = hist_dmes.GetYaxis().FindBin(0.5999)
            ptdmes_min = hist_dmes.GetXaxis().FindBin(0.5001)
            ptdmes_max = hist_dmes.GetXaxis().FindBin(23.9999)
            xsec_jpsi = hist_jpsi.Integral(yjpsi_min, yjpsi_max, ptjpsi_min, ptjpsi_max, 1, 2) * xsec / nevents
            xsec_dmes = hist_dmes.Integral(ptdmes_min, ptdmes_max, ydmes_min, ydmes_max, 1, 2) * xsec / nevents
            xsec_dmes_prompt = hist_dmes.Integral(ptdmes_min, ptdmes_max, ydmes_min, ydmes_max, 1, 1) * xsec / nevents

        df["deltaY"] = df["yJPsi"] - df["yD"]
        df["deltaPhi"] = df.apply(lambda row: abs(row.phiJPsi - row.phiD) if abs(row.phiJPsi - row.phiD) < np.pi else 2*np.pi - abs(row.phiJPsi - row.phiD), axis=1)

        df_sps = df.query(f"commonAnchestor > 0.5 and abs(pdgD) == 421{kine_cuts}")
        hist_deltay_sps_jpsiincl_dmesincl = ROOT.TH1D("hist_deltay_sps_jpsiincl_dmesincl",
                                                      ";#Delta#it{y}; d#sigma/d#Delta#it{y} (#mub)",
                                                      100, -10., 10.)
        hist_deltay_sps_jpsiincl_dmesprompt = ROOT.TH1D("hist_deltay_sps_jpsiincl_dmesprompt",
                                                        ";#Delta#it{y}; d#sigma/d#Delta#it{y} (#mub)",
                                                        100, -10., 10.)
        hist_deltaphi_sps_jpsiincl_dmesincl = ROOT.TH1D("hist_deltaphi_sps_jpsiincl_dmesincl",
                                                        ";|#Delta#varphi|; d#sigma/d|#Delta#varphi| (#mub)",
                                                        90, 0., np.pi)
        hist_deltaphi_sps_jpsiincl_dmesprompt = ROOT.TH1D("hist_deltaphi_sps_jpsiincl_dmesprompt",
                                                          ";|#Delta#varphi|; d#sigma/|d#Delta#varphi| (#mub)",
                                                          90, 0., np.pi)

        for deltay, deltaphi, dfromb in zip(df_sps["deltaY"].to_numpy(),
                                            df_sps["deltaPhi"].to_numpy(),
                                            df_sps["DfromB"].to_numpy()):
            hist_deltay_sps_jpsiincl_dmesincl.Fill(deltay)
            hist_deltaphi_sps_jpsiincl_dmesincl.Fill(deltaphi)
            if not dfromb:
                hist_deltay_sps_jpsiincl_dmesprompt.Fill(deltay)
                hist_deltaphi_sps_jpsiincl_dmesprompt.Fill(deltaphi)

        xsec_sps_jpsidmesincl = hist_deltay_sps_jpsiincl_dmesincl.Integral() * xsec / nevents
        unc_sps_incl = np.sqrt(hist_deltay_sps_jpsiincl_dmesincl.Integral()) * xsec / nevents
        xsec_sps_jpsidmesprompt = hist_deltay_sps_jpsiincl_dmesprompt.Integral() * xsec / nevents
        unc_sps_prompt = np.sqrt(hist_deltay_sps_jpsiincl_dmesprompt.Integral()) * xsec / nevents
        ratio_sps_incl = xsec_jpsi * xsec_dmes / xsec_sps_jpsidmesincl
        ratiounc_sps_incl = unc_sps_incl / xsec_sps_jpsidmesincl * ratio_sps_incl
        hist_ratio_sps_jpsiincl_dmesincl.SetBinContent(ivar+1, ratio_sps_incl)
        hist_ratio_sps_jpsiincl_dmesincl.SetBinError(ivar+1, ratiounc_sps_incl)
        ratio_sps_prompt = xsec_jpsi * xsec_dmes / xsec_sps_jpsidmesprompt
        ratiounc_sps_prompt = unc_sps_prompt / xsec_sps_jpsidmesprompt * ratio_sps_prompt
        hist_ratio_sps_jpsiincl_dmesprompt.SetBinContent(ivar+1, ratio_sps_prompt)
        hist_ratio_sps_jpsiincl_dmesprompt.SetBinError(ivar+1, ratiounc_sps_prompt)

        hist_deltay_sps_jpsiincl_dmesincl.Scale(xsec/nevents, "width")
        hist_deltay_sps_jpsiincl_dmesprompt.Scale(xsec/nevents, "width")
        hist_deltaphi_sps_jpsiincl_dmesincl.Scale(xsec/nevents, "width")
        hist_deltaphi_sps_jpsiincl_dmesprompt.Scale(xsec/nevents, "width")

        set_style(hist_deltay_sps_jpsiincl_dmesincl, COLORS[0])
        set_style(hist_deltay_sps_jpsiincl_dmesprompt, COLORS[0])
        set_style(hist_deltaphi_sps_jpsiincl_dmesincl, COLORS[0])
        set_style(hist_deltaphi_sps_jpsiincl_dmesprompt, COLORS[0])

        df_dps = df.query(f"commonAnchestor < 0.5 and abs(pdgD) == 421{kine_cuts}")
        hist_deltay_dps_jpsiincl_dmesincl = ROOT.TH1D("hist_deltay_dps_jpsiincl_dmesincl",
                                                      ";#Delta#it{y}; d#sigma/d#Delta#it{y} (#mub)",
                                                      100, -10., 10.)
        hist_deltay_dps_jpsiincl_dmesprompt = ROOT.TH1D("hist_deltay_dps_jpsiincl_dmesprompt",
                                                        ";#Delta#it{y}; d#sigma/d#Delta#it{y} (#mub)",
                                                        100, -10., 10.)
        hist_deltaphi_dps_jpsiincl_dmesincl = ROOT.TH1D("hist_deltaphi_dps_jpsiincl_dmesincl",
                                                        ";|#Delta#varphi|; d#sigma/d|#Delta#varphi| (#mub)",
                                                        90, 0., np.pi)
        hist_deltaphi_dps_jpsiincl_dmesprompt = ROOT.TH1D("hist_deltaphi_dps_jpsiincl_dmesprompt",
                                                          ";|#Delta#varphi|; d#sigma/d|#Delta#varphi| (#mub)",
                                                          90, 0., np.pi)
        for deltay, deltaphi, dfromb in zip(df_dps["deltaY"].to_numpy(),
                                            df_dps["deltaPhi"].to_numpy(),
                                            df_dps["DfromB"].to_numpy()):
            hist_deltay_dps_jpsiincl_dmesincl.Fill(deltay)
            hist_deltaphi_dps_jpsiincl_dmesincl.Fill(deltaphi)
            if not dfromb:
                hist_deltay_dps_jpsiincl_dmesprompt.Fill(deltay)
                hist_deltaphi_dps_jpsiincl_dmesprompt.Fill(deltaphi)

        xsec_dps_jpsidmesincl = hist_deltay_dps_jpsiincl_dmesincl.Integral() * xsec / nevents
        unc_dps_incl = np.sqrt(hist_deltay_dps_jpsiincl_dmesincl.Integral()) * xsec / nevents
        xsec_dps_jpsidmesprompt = hist_deltay_dps_jpsiincl_dmesprompt.Integral() * xsec / nevents
        unc_dps_prompt = np.sqrt(hist_deltay_dps_jpsiincl_dmesprompt.Integral()) * xsec / nevents
        ratio_dps_incl = xsec_jpsi * xsec_dmes / xsec_dps_jpsidmesincl
        ratio_dps_prompt = xsec_jpsi * xsec_dmes_prompt / xsec_dps_jpsidmesprompt
        ratiounc_dps_incl = unc_dps_incl / xsec_dps_jpsidmesincl * ratio_dps_incl
        ratiounc_dps_prompt = unc_dps_prompt / xsec_dps_jpsidmesprompt * ratio_dps_prompt
        unc_tot_incl = np.sqrt(hist_deltay_dps_jpsiincl_dmesincl.Integral() + xsec_sps_jpsidmesincl * nevents / xsec) * xsec / nevents
        unc_tot_prompt = np.sqrt(hist_deltay_dps_jpsiincl_dmesprompt.Integral() + xsec_sps_jpsidmesprompt * nevents / xsec) * xsec / nevents
        ratio_tot_incl = xsec_jpsi * xsec_dmes / (xsec_sps_jpsidmesincl + xsec_dps_jpsidmesincl)
        ratio_tot_prompt = xsec_jpsi * xsec_dmes_prompt / (xsec_sps_jpsidmesprompt + xsec_dps_jpsidmesprompt)
        ratiounc_tot_incl = unc_tot_incl / (xsec_sps_jpsidmesincl + xsec_dps_jpsidmesincl) * ratio_tot_incl
        ratiounc_tot_prompt = unc_tot_prompt / (xsec_sps_jpsidmesprompt + xsec_dps_jpsidmesprompt) * ratio_tot_prompt
        hist_ratio_dps_jpsiincl_dmesincl.SetBinContent(ivar+1, ratio_dps_incl)
        hist_ratio_dps_jpsiincl_dmesprompt.SetBinContent(ivar+1, ratio_dps_prompt)
        hist_ratio_tot_jpsiincl_dmesincl.SetBinContent(ivar+1, ratio_tot_incl)
        hist_ratio_tot_jpsiincl_dmesprompt.SetBinContent(ivar+1, ratio_tot_prompt)
        hist_ratio_dps_jpsiincl_dmesincl.SetBinError(ivar+1, ratiounc_dps_incl)
        hist_ratio_dps_jpsiincl_dmesprompt.SetBinError(ivar+1, ratiounc_dps_prompt)
        hist_ratio_tot_jpsiincl_dmesincl.SetBinError(ivar+1, ratiounc_tot_incl)
        hist_ratio_tot_jpsiincl_dmesprompt.SetBinError(ivar+1, ratiounc_tot_prompt)

        print(xsec_jpsi, xsec_dmes_prompt, xsec_dps_jpsidmesprompt)

        hist_deltay_dps_jpsiincl_dmesincl.Scale(xsec/nevents, "width")
        hist_deltay_dps_jpsiincl_dmesprompt.Scale(xsec/nevents, "width")
        hist_deltaphi_dps_jpsiincl_dmesincl.Scale(xsec/nevents, "width")
        hist_deltaphi_dps_jpsiincl_dmesprompt.Scale(xsec/nevents, "width")

        set_style(hist_deltay_dps_jpsiincl_dmesincl, COLORS[1])
        set_style(hist_deltay_dps_jpsiincl_dmesprompt, COLORS[1])
        set_style(hist_deltaphi_dps_jpsiincl_dmesincl, COLORS[1])
        set_style(hist_deltaphi_dps_jpsiincl_dmesprompt, COLORS[1])

        hist_deltay_tot_jpsiincl_dmesincl = hist_deltay_dps_jpsiincl_dmesincl.Clone(
            "hist_deltay_tot_jpsiincl_dmesincl")
        hist_deltay_tot_jpsiincl_dmesprompt = hist_deltay_dps_jpsiincl_dmesprompt.Clone(
            "hist_deltay_tot_jpsiincl_dmesprompt")
        hist_deltaphi_tot_jpsiincl_dmesincl = hist_deltaphi_dps_jpsiincl_dmesincl.Clone(
            "hist_deltaphi_tot_jpsiincl_dmesincl")
        hist_deltaphi_tot_jpsiincl_dmesprompt = hist_deltaphi_dps_jpsiincl_dmesprompt.Clone(
            "hist_deltaphi_tot_jpsiincl_dmesprompt")
        hist_deltay_tot_jpsiincl_dmesincl.Add(hist_deltay_sps_jpsiincl_dmesincl)
        hist_deltay_tot_jpsiincl_dmesprompt.Add(hist_deltay_sps_jpsiincl_dmesprompt)
        hist_deltaphi_tot_jpsiincl_dmesincl.Add(hist_deltaphi_sps_jpsiincl_dmesincl)
        hist_deltaphi_tot_jpsiincl_dmesprompt.Add(hist_deltaphi_sps_jpsiincl_dmesprompt)

        set_style(hist_deltay_tot_jpsiincl_dmesincl, COLORS[2])
        set_style(hist_deltay_tot_jpsiincl_dmesprompt, COLORS[2])
        set_style(hist_deltaphi_tot_jpsiincl_dmesincl, COLORS[2])
        set_style(hist_deltaphi_tot_jpsiincl_dmesprompt, COLORS[2])

        canv[variant].cd(1)
        hist_deltay_tot_jpsiincl_dmesprompt.GetYaxis().SetDecimals()
        hist_deltay_tot_jpsiincl_dmesprompt.DrawCopy()
        hist_deltay_sps_jpsiincl_dmesprompt.DrawCopy("same")
        hist_deltay_dps_jpsiincl_dmesprompt.DrawCopy("same")
        canv[variant].cd(2)
        hist_deltaphi_tot_jpsiincl_dmesprompt.GetYaxis().SetDecimals()
        hist_deltaphi_tot_jpsiincl_dmesprompt.GetYaxis().SetRangeUser(
            0., hist_deltaphi_tot_jpsiincl_dmesprompt.GetMaximum()*1.2)
        hist_deltaphi_tot_jpsiincl_dmesprompt.DrawCopy()
        hist_deltaphi_sps_jpsiincl_dmesprompt.DrawCopy("same")
        hist_deltaphi_dps_jpsiincl_dmesprompt.DrawCopy("same")

        canv[variant].SaveAs(outfile_name.replace(".root", f"_{variant}.pdf"))

        outfile.cd()
        outfile.mkdir(variant)
        outfile.cd(variant)
        canv[variant].Write()
        hist_deltay_sps_jpsiincl_dmesincl.Write()
        hist_deltay_sps_jpsiincl_dmesprompt.Write()
        hist_deltaphi_sps_jpsiincl_dmesincl.Write()
        hist_deltaphi_sps_jpsiincl_dmesprompt.Write()
        hist_deltay_dps_jpsiincl_dmesincl.Write()
        hist_deltay_dps_jpsiincl_dmesprompt.Write()
        hist_deltaphi_dps_jpsiincl_dmesincl.Write()
        hist_deltaphi_dps_jpsiincl_dmesprompt.Write()
        hist_deltay_tot_jpsiincl_dmesincl.Write()
        hist_deltay_tot_jpsiincl_dmesprompt.Write()
        hist_deltaphi_tot_jpsiincl_dmesincl.Write()
        hist_deltaphi_tot_jpsiincl_dmesprompt.Write()

    canv_ratio = ROOT.TCanvas("canv_ratio", "", 1000, 500)
    canv_ratio.SetLogy()
    hist_ratio_dps_jpsiincl_dmesprompt.GetYaxis().SetRangeUser(1.e2, 1.e5)
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
    args = parser.parse_args()

    produce_predictions(args.infile, args.outfile)
