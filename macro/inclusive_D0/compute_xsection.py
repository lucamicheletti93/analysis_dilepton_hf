"""
Compute cross section from output of cut-variation method
"""

import sys
import argparse
import numpy as np
import ROOT
from style_formatter import set_object_style


def compute_xsec(input_rawy, input_eff, input_cutvar, input_norm, had, suffix):
    """
    Main function for computation of pT-differential and pT-integrated cross sections
    """

    if had == "dplus":
        br = 0.0938
    elif had == "dstar":
        br = 0.677 * 0.03951
    elif had == "dzero":
        br = 0.03951
    else:
        print(f"ERROR: particle {had} not supported, exit")
        sys.exit()

    infile_cutvar = ROOT.TFile.Open(input_cutvar)
    hist_xsec_p = infile_cutvar.Get("hRawFracPrompt")
    hist_xsec_np = infile_cutvar.Get("hRawFracNonPrompt")
    hist_xsec_p.SetDirectory(0)
    hist_xsec_np.SetDirectory(0)
    hist_xsec_p.SetName("hist_xsec_p")
    hist_xsec_np.SetName("hist_xsec_np")
    hist_xsec_p.Sumw2(0)
    hist_xsec_np.Sumw2(0)
    set_object_style(hist_xsec_p, color=ROOT.kRed+1)
    set_object_style(hist_xsec_np, color=ROOT.kAzure+4)
    infile_cutvar.Close()

    infile_efficiencies = ROOT.TFile.Open(input_eff)
    hist_eff_prompt = infile_efficiencies.Get("hist_eff_prompt")
    hist_eff_nonprompt = infile_efficiencies.Get("hist_eff_nonprompt")
    hist_eff_prompt.SetDirectory(0)
    hist_eff_nonprompt.SetDirectory(0)
    hist_eff_prompt.Sumw2()
    hist_eff_nonprompt.Sumw2()
    infile_efficiencies.Close()

    infile_rawy = ROOT.TFile.Open(input_rawy)
    hist_rawyield = infile_rawy.Get("hist_rawyield")
    hist_rawyield.SetDirectory(0)
    hist_rawyield.Sumw2()
    infile_rawy.Close()

    infile_cutvar = ROOT.TFile.Open(input_cutvar)
    hist_xsec_p = infile_cutvar.Get("hRawFracPrompt")
    hist_xsec_np = infile_cutvar.Get("hRawFracNonPrompt")
    hist_xsec_p.SetDirectory(0)
    hist_xsec_np.SetDirectory(0)
    hist_xsec_p.SetName("hist_xsec_p")
    hist_xsec_np.SetName("hist_xsec_np")
    hist_xsec_p.Sumw2()
    hist_xsec_np.Sumw2()
    set_object_style(hist_xsec_p, color=ROOT.kRed+1)
    set_object_style(hist_xsec_np, color=ROOT.kAzure+4)
    infile_cutvar.Close()

    hist_xsec_p.Multiply(hist_rawyield)
    hist_xsec_np.Multiply(hist_rawyield)
    hist_xsec_p.Divide(hist_eff_prompt)
    hist_xsec_np.Divide(hist_eff_nonprompt)

    infile_norm = ROOT.TFile.Open(input_norm)
    hist_norm = infile_norm.Get("hist_lumi")
    lumi = hist_norm.GetBinContent(1)
    infile_cutvar.Close()

    hist_xsec_p.Scale(1./2/lumi/br, "width")
    hist_xsec_np.Scale(1./2/lumi/br, "width")
    hist_xsec_p.SetTitle(";#it{p}_{T} (GeV/#it{c});d#sigma/d#it{p}_{T} (#mub GeV^{#minus1} #it{c})")
    hist_xsec_np.SetTitle(";#it{p}_{T} (GeV/#it{c});d#sigma/d#it{p}_{T} (#mub GeV^{#minus1} #it{c})")

    nptbins = hist_xsec_p.GetNbinsX()
    pt_min = hist_xsec_p.GetBinLowEdge(1)
    pt_max = hist_xsec_p.GetXaxis().GetBinUpEdge(nptbins)
    hist_xsec_ptint_p = ROOT.TH1D(
        "hist_xsec_ptint_p", ";#it{p}_{T} (GeV/#it{c}); #sigma (#mub)", 1, pt_min, pt_max)
    hist_xsec_ptint_np = ROOT.TH1D(
        "hist_xsec_ptint_np", ";#it{p}_{T} (GeV/#it{c}); #sigma (#mub)", 1, pt_min, pt_max)
    hist_xsec_ptgrt1_p = ROOT.TH1D(
        "hist_xsec_ptgrt1_p", ";#it{p}_{T} (GeV/#it{c}); #sigma (#mub)", 1, 1., pt_max)
    hist_xsec_ptgrt1_np = ROOT.TH1D(
        "hist_xsec_ptgrt1_np", ";#it{p}_{T} (GeV/#it{c}); #sigma (#mub)", 1, 1., pt_max)
    xsec_ptint_p, xsec_ptint_np, unc_xsec_ptint_p, unc_xsec_ptint_np = 0, 0, 0, 0
    xsec_ptgrt1_p, xsec_ptgrt1_np, unc_xsec_ptgrt1_p, unc_xsec_ptgrt1_np = 0, 0, 0, 0
    for ipt in range(nptbins):
        xsec_ptint_p += hist_xsec_p.GetBinContent(ipt+1) * hist_xsec_p.GetBinWidth(ipt+1)
        xsec_ptint_np += hist_xsec_np.GetBinContent(ipt+1) * hist_xsec_p.GetBinWidth(ipt+1)
        unc_xsec_ptint_p += hist_xsec_p.GetBinError(ipt+1)**2 * hist_xsec_p.GetBinWidth(ipt+1)**2
        unc_xsec_ptint_np += hist_xsec_np.GetBinError(ipt+1)**2 * hist_xsec_p.GetBinWidth(ipt+1)**2
        if hist_xsec_p.GetBinLowEdge(ipt+1) > 0.99:
            xsec_ptgrt1_p += hist_xsec_p.GetBinContent(ipt+1) * hist_xsec_p.GetBinWidth(ipt+1)
            xsec_ptgrt1_np += hist_xsec_np.GetBinContent(ipt+1) * hist_xsec_p.GetBinWidth(ipt+1)
            unc_xsec_ptgrt1_p += hist_xsec_p.GetBinError(ipt+1)**2 * hist_xsec_p.GetBinWidth(ipt+1)**2
            unc_xsec_ptgrt1_np += hist_xsec_np.GetBinError(ipt+1)**2 * hist_xsec_p.GetBinWidth(ipt+1)**2
            
    unc_xsec_ptint_p = np.sqrt(unc_xsec_ptint_p)
    unc_xsec_ptint_np = np.sqrt(unc_xsec_ptint_np)
    hist_xsec_ptint_p.SetBinContent(1, xsec_ptint_p)
    hist_xsec_ptint_p.SetBinError(1, unc_xsec_ptint_p)
    hist_xsec_ptint_np.SetBinContent(1, xsec_ptint_np)
    hist_xsec_ptint_np.SetBinError(1, unc_xsec_ptint_np)
    unc_xsec_ptgrt1_p = np.sqrt(unc_xsec_ptgrt1_p)
    unc_xsec_ptgrt1_np = np.sqrt(unc_xsec_ptgrt1_np)
    hist_xsec_ptgrt1_p.SetBinContent(1, xsec_ptgrt1_p)
    hist_xsec_ptgrt1_p.SetBinError(1, unc_xsec_ptgrt1_p)
    hist_xsec_ptgrt1_np.SetBinContent(1, xsec_ptgrt1_np)
    hist_xsec_ptgrt1_np.SetBinError(1, unc_xsec_ptgrt1_np)

    outfile = ROOT.TFile(f"../../data_shared/{had}_xsec_pp13dot6TeV{suffix}.root", "recreate")
    hist_xsec_p.Write()
    hist_xsec_np.Write()
    hist_xsec_ptint_p.Write()
    hist_xsec_ptint_np.Write()
    hist_xsec_ptgrt1_p.Write()
    hist_xsec_ptgrt1_np.Write()
    hist_norm.Write()
    outfile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--input_rawyield", "-ir", metavar="text",
                        default="rawyields/rawyields_nocut_dzero_LHC24_JPsiD.root",
                        help="input file with raw yields")
    parser.add_argument("--input_efficiency", "-ie", metavar="text",
                        default="efficiencies/efficiencies_nocutnp_dzero_LHC24k3_trackTuner_ptSmearing1p5_JPsiD_y06.root",
                        help="input file with efficiencies")
    parser.add_argument("--input_cutvar", "-ic", metavar="text",
                        default="cutvariation/promptfrac_dzero_pp13dot6tev_LHC24_JPsiD_y06.root",
                        help="input file with cut-variation output")
    parser.add_argument("--input_normalisation", "-in", metavar="text",
                        default="../../data_shared/luminosity_dzero_LHC24_minBias_sampled.root",
                        help="input file with normalisation")
    parser.add_argument("--particle", "-p", metavar="text",
                        default="dzero",
                        help="Particle species, options [dplus, dstar, dzero]")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="_y06", help="suffix for output file")
    args = parser.parse_args()

    compute_xsec(args.input_rawyield, args.input_efficiency, args.input_cutvar,
                 args.input_normalisation, args.particle, args.suffix)
