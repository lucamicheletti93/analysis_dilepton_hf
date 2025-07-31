"""
"""

import argparse
import os
import numpy as np
import yaml
import ROOT

COLORS = [
    ROOT.kBlack,
    ROOT.kViolet+3,
    ROOT.kViolet+10,
    ROOT.kViolet+7,
    ROOT.kBlue+2,
    ROOT.kBlue,
    ROOT.kAzure+4,
    ROOT.kAzure+2,
    ROOT.kGreen+3,
    ROOT.kGreen+2,
    ROOT.kSpring-5,
    ROOT.kOrange-2,
    ROOT.kOrange+7,
    ROOT.kRed,
    ROOT.kRed+2,
    ROOT.kRed+3
]


def set_style(hist, color):
    """
    """

    hist.SetLineWidth(2)
    hist.SetMarkerStyle(ROOT.kFullCircle)
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetMarkerSize(1)


def plot_systematics(cfg_file):
    """
    """

    ROOT.gStyle.SetTitleSize(0.045, "xy")
    ROOT.gStyle.SetLabelSize(0.045, "xy")
    ROOT.gStyle.SetPadRightMargin(0.035)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadBottomMargin(0.1)
    ROOT.gStyle.SetPadTopMargin(0.035)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

    with open(cfg_file, "r") as f:
        cfg = yaml.safe_load(f)

    hist_rawy, hist_mass, hist_width, hist_effp, hist_effnp, hist_frac, hist_xsec_ptdiff, hist_xsec = ({} for _ in range(8))
    hratio_rawy,hratio_mass,hratio_width,hratio_effp,hratio_effnp,hratio_frac,hratio_xsec_ptdiff = ({} for _ in range(7))
    suffix = "_central"
    infile_eff = ROOT.TFile.Open(cfg["input"]["central"]["efficiency"])
    hist_effp[suffix] = infile_eff.Get("hist_eff_prompt")
    hist_effp[suffix].SetDirectory(0)
    hist_effp[suffix].SetName(f"{hist_effp[suffix].GetName()}{suffix}")
    hist_effnp[suffix] = infile_eff.Get("hist_eff_nonprompt")
    hist_effnp[suffix].SetDirectory(0)
    hist_effnp[suffix].SetName(f"{hist_effnp[suffix].GetName()}{suffix}")
    set_style(hist_effp[suffix], COLORS[0])
    set_style(hist_effnp[suffix], COLORS[0])
    infile_rawy = ROOT.TFile.Open(cfg["input"]["central"]["rawyield"])
    hist_rawy[suffix] = infile_rawy.Get("hist_rawyield")
    hist_rawy[suffix].SetDirectory(0)
    hist_rawy[suffix].SetName(f"{hist_rawy[suffix].GetName()}{suffix}")
    hist_width[suffix] = infile_rawy.Get("hist_sigma")
    hist_width[suffix].SetDirectory(0)
    hist_width[suffix].SetName(f"{hist_width[suffix].GetName()}{suffix}")
    hist_mass[suffix] = infile_rawy.Get("hist_mean")
    hist_mass[suffix].SetDirectory(0)
    hist_mass[suffix].SetName(f"{hist_mass[suffix].GetName()}{suffix}")
    set_style(hist_rawy[suffix], COLORS[0])
    set_style(hist_mass[suffix], COLORS[0])
    set_style(hist_width[suffix], COLORS[0])
    infile_frac = ROOT.TFile.Open(cfg["input"]["central"]["fraction"])
    hist_frac[suffix] = infile_frac.Get("hRawFracPrompt")
    hist_frac[suffix].SetDirectory(0)
    hist_frac[suffix].SetName(f"{hist_frac[suffix].GetName()}{suffix}")
    set_style(hist_frac[suffix], COLORS[0])
    infile_xsec = ROOT.TFile.Open(cfg["input"]["central"]["xsection"])
    hist_xsec_ptdiff[suffix] = infile_xsec.Get("hist_xsec_p")
    hist_xsec_ptdiff[suffix].SetDirectory(0)
    hist_xsec_ptdiff[suffix].SetName(f"{hist_xsec_ptdiff[suffix].GetName()}{suffix}")
    hist_xsec[suffix] = infile_xsec.Get("hist_xsec_ptint_p")
    hist_xsec[suffix].SetDirectory(0)
    hist_xsec[suffix].SetName(f"{hist_xsec[suffix].GetName()}{suffix}")
    set_style(hist_xsec_ptdiff[suffix], COLORS[0])

    pt_min = hist_rawy[suffix].GetBinLowEdge(1)
    pt_max = hist_rawy[suffix].GetXaxis().GetBinUpEdge(hist_rawy[suffix].GetNbinsX())
    rawy_min = hist_rawy[suffix].GetMinimum() / 10
    rawy_max = hist_rawy[suffix].GetMaximum() * 5
    xsec_cent = hist_xsec[suffix].GetBinContent(1)
    xsecs = [xsec_cent]

    input_dir = cfg["input"]["variations"]["directory"]
    suffixes = cfg["input"]["variations"]["suffixes"]
    n_files = len(suffixes)
    hist_xsec_vstrial = ROOT.TH1D("hist_xsec_vstrial", ";trial; #sigma (#mub)", n_files+1, 0.5, n_files+1.5)
    hratio_xsec_vstrial = ROOT.TH1D("hratio_xsec_vstrial", ";trial; #sigma ratio", n_files, 0.5, n_files+0.5)
    hist_xsec_vstrial.SetDirectory(0)
    hratio_xsec_vstrial.SetDirectory(0)
    set_style(hist_xsec_vstrial, COLORS[0])
    set_style(hratio_xsec_vstrial, COLORS[0])
    hist_xsec_vstrial.SetBinContent(1, xsec_cent)
    hist_xsec_vstrial.SetBinError(1, hist_xsec[suffix].GetBinError(1))
    for file_name in os.listdir(input_dir):
        if ".root" in file_name:
            if "efficiencies" in file_name:
                infile_eff = ROOT.TFile.Open(os.path.join(input_dir, file_name))
                for isuffix, suffix in enumerate(suffixes):
                    if suffix in file_name:
                        hist_effp[suffix] = infile_eff.Get("hist_eff_prompt")
                        hist_effp[suffix].SetDirectory(0)
                        hist_effp[suffix].SetName(f"{hist_effp[suffix].GetName()}{suffix}")
                        hist_effnp[suffix] = infile_eff.Get("hist_eff_nonprompt")
                        hist_effnp[suffix].SetDirectory(0)
                        hist_effnp[suffix].SetName(f"{hist_effnp[suffix].GetName()}{suffix}")
                        set_style(hist_effp[suffix], COLORS[isuffix+1])
                        set_style(hist_effnp[suffix], COLORS[isuffix+1])
                        hratio_effp[suffix] = hist_effp[suffix].Clone(
                            f"{hist_effp[suffix].GetName().replace('hist', 'hratio')}{suffix}")
                        hratio_effp[suffix].Divide(hist_effp[suffix], hist_effp["_central"], 1., 1., "B")
                        hratio_effnp[suffix] = hist_effnp[suffix].Clone(
                            f"{hist_effnp[suffix].GetName().replace('hist', 'hratio')}{suffix}")
                        hratio_effnp[suffix].Divide(hist_effnp[suffix], hist_effnp["_central"], 1., 1., "B")
                        hratio_effp[suffix].SetDirectory(0)
                        hratio_effnp[suffix].SetDirectory(0)
            elif "rawyields" in file_name:
                infile_rawy = ROOT.TFile.Open(os.path.join(input_dir, file_name))
                for isuffix, suffix in enumerate(suffixes):
                    if suffix in file_name:
                        hist_rawy[suffix] = infile_rawy.Get("hist_rawyield")
                        hist_rawy[suffix].SetDirectory(0)
                        hist_rawy[suffix].SetName(f"{hist_rawy[suffix].GetName()}{suffix}")
                        hist_width[suffix] = infile_rawy.Get("hist_sigma")
                        hist_width[suffix].SetDirectory(0)
                        hist_width[suffix].SetName(f"{hist_width[suffix].GetName()}{suffix}")
                        hist_mass[suffix] = infile_rawy.Get("hist_mean")
                        hist_mass[suffix].SetDirectory(0)
                        hist_mass[suffix].SetName(f"{hist_mass[suffix].GetName()}{suffix}")
                        set_style(hist_rawy[suffix], COLORS[isuffix+1])
                        set_style(hist_mass[suffix], COLORS[isuffix+1])
                        set_style(hist_width[suffix], COLORS[isuffix+1])
                        hratio_rawy[suffix] = hist_rawy[suffix].Clone(
                            f"{hist_rawy[suffix].GetName().replace('hist', 'hratio')}{suffix}")
                        hratio_rawy[suffix].Divide(hist_rawy[suffix], hist_rawy["_central"], 1., 1., "B")
                        hratio_mass[suffix] = hist_mass[suffix].Clone(
                            f"{hist_mass[suffix].GetName().replace('hist', 'hratio')}{suffix}")
                        hratio_mass[suffix].Divide(hist_mass[suffix], hist_mass["_central"], 1., 1., "B")
                        hratio_width[suffix] = hist_width[suffix].Clone(
                            f"{hist_width[suffix].GetName().replace('hist', 'hratio')}{suffix}")
                        hratio_width[suffix].Divide(hist_width[suffix], hist_width["_central"], 1., 1., "B")
                        hratio_rawy[suffix].SetDirectory(0)
                        hratio_mass[suffix].SetDirectory(0)
                        hratio_width[suffix].SetDirectory(0)
            elif "promptfrac" in file_name:
                infile_frac = ROOT.TFile.Open(os.path.join(input_dir, file_name))
                for isuffix, suffix in enumerate(suffixes):
                    if suffix in file_name:
                        hist_frac[suffix] = infile_frac.Get("hRawFracPrompt")
                        hist_frac[suffix].SetDirectory(0)
                        hist_frac[suffix].SetName(f"{hist_frac[suffix].GetName()}{suffix}")
                        set_style(hist_frac[suffix], COLORS[isuffix+1])
                        hratio_frac[suffix] = hist_frac[suffix].Clone(
                            f"{hist_frac[suffix].GetName().replace('h', 'hRatio')}{suffix}")
                        hratio_frac[suffix].Divide(hist_frac[suffix], hist_frac["_central"], 1., 1., "B")
                        hratio_frac[suffix].SetDirectory(0)
            elif "dzero_xsec" in file_name:
                infile_xsec = ROOT.TFile.Open(os.path.join(input_dir, file_name))
                for isuffix, suffix in enumerate(suffixes):
                    if suffix in file_name:
                        hist_xsec_ptdiff[suffix] = infile_xsec.Get("hist_xsec_p")
                        hist_xsec_ptdiff[suffix].SetDirectory(0)
                        hist_xsec_ptdiff[suffix].SetName(f"{hist_xsec_ptdiff[suffix].GetName()}{suffix}")
                        hist_xsec[suffix] = infile_xsec.Get("hist_xsec_ptint_p")
                        hist_xsec[suffix].SetDirectory(0)
                        hist_xsec[suffix].SetName(f"{hist_xsec[suffix].GetName()}{suffix}")
                        set_style(hist_xsec_ptdiff[suffix], COLORS[isuffix+1])
                        hist_xsec_vstrial.SetBinContent(isuffix+2, hist_xsec[suffix].GetBinContent(1))
                        hist_xsec_vstrial.SetBinError(isuffix+2, hist_xsec[suffix].GetBinError(1))
                        hratio_xsec_ptdiff[suffix] = hist_xsec_ptdiff[suffix].Clone(
                            f"{hist_xsec_ptdiff[suffix].GetName().replace('h', 'hRatio')}{suffix}")
                        hratio_xsec_ptdiff[suffix].Divide(hist_xsec_ptdiff[suffix], hist_xsec_ptdiff["_central"], 1., 1., "B")
                        hratio_xsec_ptdiff[suffix].SetDirectory(0)
                        xsecs.append(hist_xsec[suffix].GetBinContent(1))
                        hratio_xsec_vstrial.SetBinContent(isuffix+1, hist_xsec[suffix].GetBinContent(1)/xsec_cent)
                        hratio_xsec_vstrial.SetBinError(isuffix+1, 1.e-20)

    # rms = np.std(xsecs)

    canv_syst = ROOT.TCanvas("canv_syst", "", 1920, 1080)
    leg = ROOT.TLegend(0.2, 0.16, 1.0, 0.8)
    leg.SetTextSize(0.035)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetNColumns(cfg["input"]["variations"]["n_cols_leg"])
    canv_syst.Divide(4, 2)
    canv_syst.cd(1).DrawFrame(pt_min, rawy_min, pt_max, rawy_max, ";#it{p}_{T} (GeV/#it{c}); raw yield")
    canv_syst.cd(1).SetLogy()
    leg.AddEntry(hist_rawy["_central"], "default", "lp")
    for isuffix, suffix in enumerate(suffixes):
        hist_rawy[suffix].DrawCopy("esame")
        leg.AddEntry(hist_rawy[suffix], cfg["input"]["variations"]["labels"][isuffix], "lp")
    canv_syst.cd(2).DrawFrame(pt_min, 1.85, pt_max, 1.87, ";#it{p}_{T} (GeV/#it{c}); peak mean (GeV/#it{c}^{2})")
    for suffix in hist_mass:
        hist_mass[suffix].DrawCopy("esame")
    canv_syst.cd(3).DrawFrame(pt_min, 0., pt_max, 0.075, ";#it{p}_{T} (GeV/#it{c});  peak width (GeV/#it{c}^{2})")
    for suffix in hist_width:
        hist_width[suffix].DrawCopy("esame")
    canv_syst.cd(4).DrawFrame(pt_min, 1.e-3, pt_max, 1., ";#it{p}_{T} (GeV/#it{c}); prompt efficiency")
    canv_syst.cd(4).SetLogy()
    for suffix in hist_effp:
        hist_effp[suffix].DrawCopy("esame")
    canv_syst.cd(5).DrawFrame(pt_min, 1.e-3, pt_max, 1., ";#it{p}_{T} (GeV/#it{c}); non-prompt efficiency")
    canv_syst.cd(5).SetLogy()
    for suffix in hist_effnp:
        hist_effnp[suffix].DrawCopy("esame")
    canv_syst.cd(6).DrawFrame(pt_min, 0., pt_max, 1., ";#it{p}_{T} (GeV/#it{c}); prompt fraction")
    for suffix in hist_frac:
        hist_frac[suffix].DrawCopy("esame")
    leg.Draw()
    canv_syst.cd(7).DrawFrame(pt_min, 1.e-2, pt_max, 1.e3, ";#it{p}_{T} (GeV/#it{c}); (d#sigma/d#it{p}_{T})_{|#it{y}|<0.6} (#mub #it{c} GeV^{#minus1})")
    canv_syst.cd(7).SetLogy()
    for suffix in hist_xsec_ptdiff:
        hist_xsec_ptdiff[suffix].DrawCopy("esame")
    canv_syst.cd(8)
    hist_xsec_vstrial.DrawCopy()
    canv_syst.SaveAs(os.path.join(cfg["output"]["directory"], cfg["output"]["file_name"]))

    canv_ratios = ROOT.TCanvas("canv_ratios", "", 1920, 1080)
    canv_ratios.Divide(4, 2)
    canv_ratios.cd(1).DrawFrame(pt_min, 0., pt_max, 2., ";#it{p}_{T} (GeV/#it{c}); raw yield ratio")
    for suffix in hratio_rawy:
        hratio_rawy[suffix].DrawCopy("esame")
    canv_ratios.cd(2).DrawFrame(pt_min, 0.99, pt_max, 1.01, ";#it{p}_{T} (GeV/#it{c}); peak mean ratio")
    for suffix in hratio_mass:
        hratio_mass[suffix].DrawCopy("esame")
    canv_ratios.cd(3).DrawFrame(pt_min, 0.7, pt_max, 1.3, ";#it{p}_{T} (GeV/#it{c});  peak width ratio")
    for suffix in hratio_width:
        hratio_width[suffix].DrawCopy("esame")
    canv_ratios.cd(4).DrawFrame(pt_min, 0., pt_max, 2., ";#it{p}_{T} (GeV/#it{c}); prompt efficiency ratio")
    for suffix in hratio_effp:
        hratio_effp[suffix].DrawCopy("esame")
    canv_ratios.cd(5).DrawFrame(pt_min, 0., pt_max, 2., ";#it{p}_{T} (GeV/#it{c}); non-prompt efficiency ratio")
    for suffix in hratio_effnp:
        hratio_effnp[suffix].DrawCopy("esame")
    canv_ratios.cd(6).DrawFrame(pt_min, 0.8, pt_max, 1.2, ";#it{p}_{T} (GeV/#it{c}); prompt fraction ratio")
    for suffix in hratio_frac:
        hratio_frac[suffix].DrawCopy("esame")
    canv_ratios.cd(7).DrawFrame(pt_min, 0.7, pt_max, 1.3, ";#it{p}_{T} (GeV/#it{c}); (d#sigma/d#it{p}_{T})_{|#it{y}|<0.6} ratio")
    for suffix in hratio_xsec_ptdiff:
        hratio_xsec_ptdiff[suffix].DrawCopy("esame")
    canv_ratios.cd(8)
    hratio_xsec_vstrial.DrawCopy()
    canv_ratios.SaveAs(os.path.join(cfg["output"]["directory"], cfg["output"]["file_name"].replace(".pdf", "_ratios.pdf")))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("-c", "--cfg_file", metavar="text", required=True,
                        help="yaml config file")
    args = parser.parse_args()

    plot_systematics(args.cfg_file)
