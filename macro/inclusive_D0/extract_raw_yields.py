"""
Script to perform signal extraction as a function of pT for several BDT output scores for cut-variation method
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import yaml
import pdg
import uproot
import ROOT
from fit_utils import perform_fit

def check_config_consistency(cfg):
    """
    
    """
    pt_mins = cfg["pt_mins"]
    pt_maxs = cfg["pt_maxs"]
    mass_mins = cfg["fit"]["mass_mins"]
    mass_maxs = cfg["fit"]["mass_maxs"]
    sgn_funcs = cfg["fit"]["sgn_funcs"]
    bkg_funcs = cfg["fit"]["bkg_funcs"]
    bdt_bkg_cuts = cfg["bdt_cuts"]["bkg"]

    n = len(pt_mins)

    if len(pt_maxs) != n:
        print("pt_mins and pt_maxs have different length, check your config!")
        sys.exit()
    if len(mass_mins) + len(mass_maxs) + len(sgn_funcs) + len(bkg_funcs) != 4*n:
        print("fit configuration length doesn't match number of pt bins, check your config!")
        sys.exit()
    if isinstance(bdt_bkg_cuts, list) and len(bdt_bkg_cuts) != n:
        print("number of bkg cuts must be a single value or have the same "
              "length of the number of pt bins, check your config!")
        sys.exit()

# function to perform fits
def fit(input_config):
    """
    Method for fitting
    """

    pdg_api = pdg.connect()

    with open(input_config, "r") as f:
        cfg = yaml.safe_load(f)

    if cfg["hadron"] not in ["dstar", "dplus", "dzero"]:
        print(f"ERROR: {cfg['hadron']} not supported, exit")
        sys.exit()

    check_config_consistency(cfg)
    outputdir = cfg["output"]["rawyields"]["directory"]
    suffix = cfg["output"]["rawyields"]["suffix"]

    pt_mins = cfg["pt_mins"]
    pt_maxs = cfg["pt_maxs"]
    pt_limits = pt_mins.copy()
    pt_limits.append(pt_maxs[-1])
    pt_limits = np.array(pt_limits, np.float64)

    infile_name = os.path.join(outputdir, f"hist_mass{suffix}.root")
    infile_name_corrbkg = os.path.join(outputdir, f"hist_bkgtempl{suffix}.root")

    mass_mins = cfg["fit"]["mass_mins"]
    mass_maxs = cfg["fit"]["mass_maxs"]
    sgn_funcs = cfg["fit"]["sgn_funcs"]
    bkg_funcs = cfg["fit"]["bkg_funcs"]
    fix_frac_bkgcorr = cfg["fit"]["fix_frac_bkgcorr"]
    fix_tailpars_from_mc = cfg["fit"]["fix_tailpars_from_mc"]
    bdt_np_mins = cfg["bdt_cuts"]["nonprompt"]

    # we first fit the high significance cases
    outfile_name_nocut = os.path.join(outputdir, f"rawyields_nocut{suffix}.root")
    outfile_nocut = ROOT.TFile(outfile_name_nocut, "recreate")
    outfile_nocut.Close()

    hist_rawyield_nocut = ROOT.TH1D(
        "hist_rawyield", ";#it{p}_{T} (GeV/#it{c}); raw yield", len(pt_mins), pt_limits)
    hist_sigma_nocut = ROOT.TH1D(
        "hist_sigma", ";#it{p}_{T} (GeV/#it{c}); #sigma (GeV/#it{c}^{2})", len(pt_mins), pt_limits)
    hist_mean_nocut = ROOT.TH1D(
        "hist_mean", ";#it{p}_{T} (GeV/#it{c}); #mu (GeV/#it{c}^{2})", len(pt_mins), pt_limits)
    hist_alphal_nocut = ROOT.TH1D(
        "hist_alphal", ";#it{p}_{T} (GeV/#it{c}); #alpha_{l}", len(pt_mins), pt_limits)
    hist_alphar_nocut = ROOT.TH1D(
        "hist_alphar", ";#it{p}_{T} (GeV/#it{c}); #alpha_{r}", len(pt_mins), pt_limits)
    hist_nl_nocut = ROOT.TH1D(
        "hist_nl", ";#it{p}_{T} (GeV/#it{c}); #it{n}_{l}", len(pt_mins), pt_limits)
    hist_nr_nocut = ROOT.TH1D(
        "hist_nr", ";#it{p}_{T} (GeV/#it{c}); #it{n}_{r}", len(pt_mins), pt_limits)

    hist_rawyield_nocut.SetDirectory(0)
    hist_sigma_nocut.SetDirectory(0)
    hist_mean_nocut.SetDirectory(0)
    hist_alphal_nocut.SetDirectory(0)
    hist_alphar_nocut.SetDirectory(0)
    hist_nl_nocut.SetDirectory(0)
    hist_nr_nocut.SetDirectory(0)

    hist_rawyield_cutvar = []
    for icut, _ in enumerate(bdt_np_mins):
        hist_rawyield_cutvar.append(
            ROOT.TH1D("hist_rawyield", ";#it{p}_{T} (GeV/#it{c}); raw yield", len(pt_mins), pt_limits))
        hist_rawyield_cutvar[icut].SetDirectory(0)

    fitter = []
    for ipt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
        mass = 0.
        name_bkg_template = ""
        if cfg["hadron"] == "dstar":
            mass = pdg_api.get_particle_by_mcid(413).mass - pdg_api.get_particle_by_mcid(421).mass
        elif cfg["hadron"] == "dplus":
            mass = pdg_api.get_particle_by_mcid(411).mass
            name_bkg_template = f"hist_templ_pt{pt_min:.1f}_{pt_max:.1f}_smooth"
        elif cfg["hadron"] == "dzero":
            mass = pdg_api.get_particle_by_mcid(421).mass
            name_bkg_template = f"hist_templ_refl_pt{pt_min:.1f}_{pt_max:.1f}"
        fitter_data, fitter_mc = perform_fit(infile_name,
                                             f"hist_mass_pt{pt_min:.1f}_{pt_max:.1f}_nocutnp",
                                             f"{cfg['hadron']}_pt{pt_min:.1f}_{pt_max:.1f}_nocutnp",
                                             sgn_funcs[ipt],
                                             bkg_funcs[ipt],
                                             mass_mins[ipt],
                                             mass_maxs[ipt],
                                             mass,
                                             cfg["hadron"],
                                             infile_name_corrbkg,
                                             name_bkg_template,
                                             f"hist_signal_pt{pt_min:.1f}_{pt_max:.1f}",
                                             fix_frac_bkgcorr[ipt],
                                             fix_tailpars_from_mc[ipt])
        fitter.append(fitter_data)

        pars_tofix = {0 : {}}
        if fitter[ipt].get_fit_result.converged:
            rawyield = fitter[ipt].get_raw_yield(0)
            sigma = fitter[ipt].get_signal_parameter(0, "sigma")
            mean = fitter[ipt].get_mass(0)
            pars_tofix[0]["sigma"] = sigma[0]
            pars_tofix[0]["mu"] = mean[0]
            hist_rawyield_nocut.SetBinContent(ipt+1, rawyield[0])
            hist_rawyield_nocut.SetBinError(ipt+1, rawyield[1])
            hist_sigma_nocut.SetBinContent(ipt+1, sigma[0])
            hist_sigma_nocut.SetBinError(ipt+1, sigma[1])
            hist_mean_nocut.SetBinContent(ipt+1, mean[0])
            hist_mean_nocut.SetBinError(ipt+1, mean[1])
            if sgn_funcs[ipt] == "doublecb":
                alphal = fitter[ipt].get_signal_parameter(0, "alphal")
                alphar = fitter[ipt].get_signal_parameter(0, "alphar")
                nl = fitter[ipt].get_signal_parameter(0, "nl")
                nr = fitter[ipt].get_signal_parameter(0, "nr")
                pars_tofix[0]["alphal"] = alphal[0]
                pars_tofix[0]["alphar"] = alphar[0]
                pars_tofix[0]["nl"] = nl[0]
                pars_tofix[0]["nr"] = nr[0]
                hist_alphal_nocut.SetBinContent(ipt+1, alphal[0])
                hist_alphal_nocut.SetBinError(ipt+1, alphal[1])
                hist_alphar_nocut.SetBinContent(ipt+1, alphar[0])
                hist_alphar_nocut.SetBinError(ipt+1, alphar[1])
                hist_nl_nocut.SetBinContent(ipt+1, nl[0])
                hist_nl_nocut.SetBinError(ipt+1, nl[1])
                hist_nr_nocut.SetBinContent(ipt+1, nr[0])
                hist_nr_nocut.SetBinError(ipt+1, nr[1])
            elif sgn_funcs[ipt] == "doublecbsymm":
                alpha = fitter[ipt].get_signal_parameter(0, "alpha")
                n = fitter[ipt].get_signal_parameter(0, "n")
                pars_tofix[0]["alpha"] = alpha[0]
                pars_tofix[0]["n"] = n[0]
                hist_alphal_nocut.SetBinContent(ipt+1, alpha[0])
                hist_alphal_nocut.SetBinError(ipt+1, alpha[1])
                hist_alphar_nocut.SetBinContent(ipt+1, alpha[0])
                hist_alphar_nocut.SetBinError(ipt+1, alpha[1])
                hist_nl_nocut.SetBinContent(ipt+1, n[0])
                hist_nl_nocut.SetBinError(ipt+1, n[1])
                hist_nr_nocut.SetBinContent(ipt+1, n[0])
                hist_nr_nocut.SetBinError(ipt+1, n[1])
            if mass_maxs[ipt] > 1.99 and cfg["hadron"] == "dplus":
                pars_tofix[1] = {}
                pars_tofix[1]["sigma"] = fitter[ipt].get_signal_parameter(1, "sigma")[0]
                pars_tofix[1]["mu"] = fitter[ipt].get_signal_parameter(1, "mu")[0]

            fitter[ipt].dump_to_root(outfile_name_nocut,
                                     option="update",
                                     suffix=f"_pt{pt_min:.1f}_{pt_max:.1f}_nocutnp")
            xaxis = ""
            if cfg["hadron"] == "dstar":
                xaxis = r"$M(\mathrm{K\pi\pi}) - M(\mathrm{K\pi})$ (GeV/$c^{2}$)"
            elif cfg["hadron"] == "dplus":
                xaxis = r"$M(\mathrm{K\pi\pi})$ (GeV/$c^{2}$)"
            if fitter_mc is not None:
                fig_mc, _ = fitter_mc.plot_mass_fit(style="ATLAS", figsize=(8, 8), axis_title=xaxis)
                fig_mc.savefig(
                    os.path.join(outputdir, f"massfit{suffix}_mc_pt{pt_min:.1f}_{pt_max:.1f}_nocutnp.pdf")
                )

            fig, _ = fitter[ipt].plot_mass_fit(style="ATLAS", figsize=(8, 8), axis_title=xaxis)
            figres = fitter[ipt].plot_raw_residuals(figsize=(8, 8), style="ATLAS")

            fig.savefig(
                os.path.join(outputdir, f"massfit{suffix}_pt{pt_min:.1f}_{pt_max:.1f}_nocutnp.pdf")
            )
            figres.savefig(
                os.path.join(outputdir, f"massfitres{suffix}_pt{pt_min:.1f}_{pt_max:.1f}_nocutnp.pdf")
            )

        fitter_pt_cutvar = []
        for icut, bdt_np_min in enumerate(bdt_np_mins):
            fitter_pt_cutvar.append(
                perform_fit(infile_name,
                            f"hist_mass_pt{pt_min:.1f}_{pt_max:.1f}_bdtnp{bdt_np_min:.2f}",
                            f"{cfg['hadron']}_pt{pt_min:.1f}_{pt_max:.1f}_bdtnp{bdt_np_min:.2f}",
                            sgn_funcs[ipt],
                            bkg_funcs[ipt],
                            mass_mins[ipt],
                            mass_maxs[ipt],
                            mass,
                            cfg["hadron"],
                            infile_name_corrbkg,
                            name_bkg_template,
                            f"hist_signal_pt{pt_min:.1f}_{pt_max:.1f}",
                            fix_frac_bkgcorr[ipt],
                            False,
                            pars_tofix)
            )
            if fitter_pt_cutvar[icut][0].get_fit_result.converged:
                if cfg["hadron"] == "dstar":
                    rawyield = fitter_pt_cutvar[icut][0].get_raw_yield_bincounting(0, min=0.14, max=0.16)
                else:
                    rawyield = fitter_pt_cutvar[icut][0].get_raw_yield_bincounting(0) # 3sigma by default
                hist_rawyield_cutvar[icut].SetBinContent(ipt+1, rawyield[0])
                hist_rawyield_cutvar[icut].SetBinError(ipt+1, rawyield[1])

    for icut, bdt_np_min in enumerate(bdt_np_mins):
        outfile_name_cutvar = os.path.join(
            outputdir, f"rawyields_bdtnp{bdt_np_min:0.2f}{suffix}.root")
        outfile_cutvar = ROOT.TFile(outfile_name_cutvar, "recreate")
        hist_rawyield_cutvar[icut].Write()
        outfile_cutvar.Close()

    outfile_nocut = ROOT.TFile(outfile_name_nocut, "update")
    hist_rawyield_nocut.Write()
    hist_sigma_nocut.Write()
    hist_mean_nocut.Write()
    hist_alphal_nocut.Write()
    hist_alphar_nocut.Write()
    hist_nl_nocut.Write()
    hist_nr_nocut.Write()
    outfile_nocut.Close()


# function to project the sparse
def project(input_config):
    """

    """

    with open(input_config, "r") as f:
        cfg = yaml.safe_load(f)

    if cfg["hadron"] not in ["dstar", "dplus", "dzero"]:
        print(f"ERROR: {cfg['hadron']} not supported, exit")
        sys.exit()

    infile_name = cfg["input"]["data"]
    outputdir = cfg["output"]["rawyields"]["directory"]
    suffix = cfg["output"]["rawyields"]["suffix"]

    mass_axis = 0
    pt_axis = 1
    bdtbkg_axis = 2
    bdtnp_axis = 3
    ncls_its_axis = None
    ncls_tpc_axis = None
    infile = ROOT.TFile.Open(infile_name)
    if cfg["hadron"] != "dzero":
        sparse = infile.Get("hData")
    else:
        sparse = infile.Get("hf-task-d0/hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type")
        mass_axis = 3
        pt_axis = 4
        bdtbkg_axis = 0
        bdtnp_axis = 2
        refl_axis = 6
        ncls_its_axis = 7
        ncls_tpc_axis = 8
    infile.Close()

    pt_mins = cfg["pt_mins"]
    pt_maxs = cfg["pt_maxs"]

    bdt_bkg_cuts = cfg["bdt_cuts"]["bkg"]
    if not isinstance(bdt_bkg_cuts, list):
        bdt_bkg_cuts = [bdt_bkg_cuts]*len(pt_mins)

    outfile = ROOT.TFile(os.path.join(outputdir, f"hist_mass{suffix}.root"), "recreate")

    if cfg["trk_cuts"]["apply"]:
        if cfg["hadron"] != "dzero":
            print("ERROR: Track cuts not implemented! Exit")
            sys.exit()
        ncls_its_bin_min = sparse.GetAxis(ncls_its_axis).FindBin(cfg["trk_cuts"]["ncls_its"]*1.001)
        ncls_tpc_bin_min = sparse.GetAxis(ncls_tpc_axis).FindBin(cfg["trk_cuts"]["ncls_tpc"]*1.001)
        sparse.GetAxis(ncls_its_axis).SetRange(ncls_its_bin_min, 100000)
        sparse.GetAxis(ncls_tpc_axis).SetRange(ncls_tpc_bin_min, 100000)

    histos_pt, histos_pt_cutvar = [], []
    for ipt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
        print(f'projecting in pt range {pt_min, pt_max}')
        bdt_bkg_bin_max = sparse.GetAxis(bdtbkg_axis).FindBin(bdt_bkg_cuts[ipt]*0.999)
        sparse.GetAxis(bdtbkg_axis).SetRange(1, bdt_bkg_bin_max)
        pt_bin_min = sparse.GetAxis(pt_axis).FindBin(pt_min*1.001)
        pt_bin_max = sparse.GetAxis(pt_axis).FindBin(pt_max*0.999)
        sparse.GetAxis(pt_axis).SetRange(pt_bin_min, pt_bin_max)
        if cfg["hadron"] == "dzero":
            sparse.GetAxis(refl_axis).SetRange(1, 2)
        histos_pt.append(sparse.Projection(mass_axis))
        histos_pt[ipt].SetName(f"hist_mass_pt{pt_min:.1f}_{pt_max:.1f}_nocutnp")
        outfile.cd()
        histos_pt[ipt].Write()
        histos_pt_cutvar.append([])
        for bdt_np_min in cfg["bdt_cuts"]["nonprompt"]:
            bdt_bkg_np_min = sparse.GetAxis(bdtnp_axis).FindBin(bdt_np_min*1.001)
            sparse.GetAxis(bdtnp_axis).SetRange(bdt_bkg_np_min, sparse.GetAxis(bdtnp_axis).GetNbins()+1)
            histos_pt_cutvar[ipt].append(sparse.Projection(mass_axis))
            histos_pt_cutvar[ipt][-1].SetName(f"hist_mass_pt{pt_min:.1f}_{pt_max:.1f}_bdtnp{bdt_np_min:.2f}")
            outfile.cd()
            histos_pt_cutvar[ipt][-1].Write()
            sparse.GetAxis(bdtnp_axis).SetRange(-1, -1)
    outfile.Close()


def get_templates(input_config):
    """
    """

    with open(input_config, "r") as f:
        cfg = yaml.safe_load(f)

    if cfg["hadron"] not in ["dplus", "dzero"]:
        print(f"ERROR: {cfg['hadron']} not supported for bkg templates, exit")
        sys.exit()

    pt_mins = cfg["pt_mins"]
    pt_maxs = cfg["pt_maxs"]
    pt_lims = pt_mins.copy()
    pt_lims.append(pt_maxs[-1])
    infile_names = cfg["input"]["corrbkg"]
    outputdir = cfg["output"]["rawyields"]["directory"]
    suffix = cfg["output"]["rawyields"]["suffix"]

    if infile_names is None:
        print("No templates to project, continue")
        return

    if cfg["hadron"] == "dplus":

        if cfg["trk_cuts"]["apply"]:
            print("ERROR: Track cuts not implemented! Exit")
            sys.exit()

        df = pd.DataFrame()
        for infile_name in infile_names:
            infile = ROOT.TFile.Open(infile_name)
            for key in infile.GetListOfKeys():
                df = pd.concat(
                    [df, uproot.open(infile_name)[f"{key.GetName()}/O2hfcanddplite"].arrays(library="pd")])

        df_bkg = df.query("abs(fFlagMcMatchRec) == 4")
        df_signal = df.query("abs(fFlagMcMatchRec) == 1")

        hist_templ_dplustokkpi, hist_templ_dstokkpi, hist_templ, hist_signal = [], [], [], []
        hist_frac_bkg_to_signal = ROOT.TH1D("hist_frac_bkg_to_signal",
                                            ";#it{p}_{T} (GeV/#it{c});bkg corr / signal",
                                            len(pt_mins), np.array(pt_lims, dtype=np.float64))
        for ipt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
            hist_templ.append(ROOT.TH1D(
                f"hist_templ_pt{pt_min:.1f}_{pt_max:.1f}",
                ";#it{M}(K#pi#pi) (GeV/#it{c})", 600, 1.67, 2.27))
            hist_templ_dplustokkpi.append(ROOT.TH1D(
                f"hist_templ_dplustokkpi_pt{pt_min:.1f}_{pt_max:.1f}",
                ";#it{M}(K#pi#pi) (GeV/#it{c})", 600, 1.67, 2.27))
            hist_templ_dstokkpi.append(ROOT.TH1D(
                f"hist_templ_dstokkpi_pt{pt_min:.1f}_{pt_max:.1f}",
                ";#it{M}(K#pi#pi) (GeV/#it{c})", 600, 1.67, 2.27))
            hist_signal.append(ROOT.TH1D(
                f"hist_signal_pt{pt_min:.1f}_{pt_max:.1f}",
                ";#it{M}(K#pi#pi) (GeV/#it{c})", 600, 1.67, 2.27))

            df_pt_dplustokkpi = df_bkg.query(f"{pt_min} < fPt < {pt_max} and abs(fFlagMcDecayChanRec) > 2")
            for mass in df_pt_dplustokkpi["fM"].to_numpy():
                hist_templ_dplustokkpi[ipt].Fill(mass-0.005)
            df_pt_dstokkpi = df_bkg.query(f"{pt_min} < fPt < {pt_max} and abs(fFlagMcDecayChanRec) <= 2")
            for mass in df_pt_dstokkpi["fM"].to_numpy():
                hist_templ_dstokkpi[ipt].Fill(mass-0.005)
            df_pt_signal = df_signal.query(f"{pt_min} < fPt < {pt_max}")
            for mass in df_pt_signal["fM"].to_numpy():
                hist_signal[ipt].Fill(mass-0.005)

            hist_templ[ipt].Add(hist_templ_dplustokkpi[ipt],
                                9.68e-3 / 0.0752 * (0.0752 + 0.0156 + 0.0104 + 0.0752))
            # Ds/D+ is underestimated in pythia CRMode2
            hist_templ[ipt].Add(hist_templ_dstokkpi[ipt], 5.37e-2 * 1.25)
            hist_signal[ipt].Scale(9.38e-2 / (0.0752 + 0.0156 + 0.0104) * (0.0752 + 0.0156 + 0.0104 + 0.0752))

            hist_frac_bkg_to_signal.SetBinContent(ipt+1,
                                                hist_templ[ipt].Integral() / hist_signal[ipt].Integral())

        outfile = ROOT.TFile(os.path.join(outputdir, f"hist_bkgtempl{suffix}.root"), "recreate")
        hist_frac_bkg_to_signal.Write()
        for hist in hist_templ_dplustokkpi:
            hist.Write()
        for hist in hist_templ_dstokkpi:
            hist.Write()
        for hist in hist_signal:
            hist.Write()
        for hist in hist_templ:
            hist.Write()
            hist_smooth = hist.Clone(f"{hist.GetName()}_smooth")
            hist_smooth.Smooth(100)
            hist_smooth.Write()
        outfile.Close()
    else:
        bdt_bkg_cuts = cfg["bdt_cuts"]["bkg"]
        if not isinstance(bdt_bkg_cuts, list):
            bdt_bkg_cuts = [bdt_bkg_cuts]*len(pt_mins)

        sparse_refl, sparse_signal = None, None
        hist_frac_bkg_to_signal = ROOT.TH1D("hist_frac_bkg_to_signal",
                                            ";#it{p}_{T} (GeV/#it{c});bkg corr / signal",
                                            len(pt_mins), np.array(pt_lims, dtype=np.float64))
        for ifile, infile_name in enumerate(infile_names):
            infile = ROOT.TFile.Open(infile_name)
            sparse_recoall = infile.Get("hf-task-d0/hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type")
            sparse_recoall.GetAxis(6).SetRange(3, 4) # only reflected signal
            if cfg["trk_cuts"]["apply"]:
                if cfg["hadron"] != "dzero":
                    print("ERROR: Track cuts not implemented! Exit")
                    sys.exit()
                ncls_its_bin_min = sparse_recoall.GetAxis(10).FindBin(cfg["trk_cuts"]["ncls_its"]*1.001)
                ncls_tpc_bin_min = sparse_recoall.GetAxis(11).FindBin(cfg["trk_cuts"]["ncls_tpc"]*1.001)
                sparse_recoall.GetAxis(10).SetRange(ncls_its_bin_min, 100000)
                sparse_recoall.GetAxis(11).SetRange(ncls_tpc_bin_min, 100000)

            if ifile == 0:
                sparse_refl = sparse_recoall.Projection(5, np.array([0, 1, 2, 3, 4], dtype=np.int32))
            else:
                sparse_refl.Add(sparse_recoall.Projection(5, np.array([0, 1, 2, 3, 4], dtype=np.int32)))
            sparse_recoall.GetAxis(6).SetRange(1, 2) # only signal
            if ifile == 0:
                sparse_signal = sparse_recoall.Projection(5, np.array([0, 1, 2, 3, 4], dtype=np.int32))
            else:
                sparse_signal.Add(sparse_recoall.Projection(5, np.array([0, 1, 2, 3, 4], dtype=np.int32)))

        hist_templ_refl, hist_signal = [], []
        for ipt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
            pt_bin_min = sparse_refl.GetAxis(4).FindBin(pt_min*1.0001)
            pt_bin_max = sparse_refl.GetAxis(4).FindBin(pt_max*0.9999)
            sparse_refl.GetAxis(4).SetRange(pt_bin_min, pt_bin_max)
            sparse_signal.GetAxis(4).SetRange(pt_bin_min, pt_bin_max)
            bdt_bin_max = sparse_refl.GetAxis(0).FindBin(bdt_bkg_cuts[ipt]*0.9999)
            sparse_refl.GetAxis(0).SetRange(1, bdt_bin_max)
            sparse_signal.GetAxis(0).SetRange(1, bdt_bin_max)
            hist_templ_refl.append(sparse_refl.Projection(3))
            hist_signal.append(sparse_signal.Projection(3))
            hist_templ_refl[ipt].SetName(f"hist_templ_refl_pt{pt_min:.1f}_{pt_max:.1f}")
            hist_signal[ipt].SetName(f"hist_signal_pt{pt_min:.1f}_{pt_max:.1f}")
            hist_frac_bkg_to_signal.SetBinContent(ipt+1, hist_templ_refl[ipt].Integral()/hist_signal[ipt].Integral())

        outfile = ROOT.TFile(os.path.join(outputdir, f"hist_bkgtempl{suffix}.root"), "recreate")
        hist_frac_bkg_to_signal.Write()
        for hist in hist_templ_refl:
            hist.Write()
        for hist in hist_signal:
            hist.Write()
        outfile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--cfg_file", "-c", metavar="text",
                        default="config.yml", help="config file")
    parser.add_argument("--project", "-p", action="store_true",
                        default=False, help="enable projection")
    parser.add_argument("--fit", "-f", action="store_true",
                        default=False, help="enable fit w/o cut")
    parser.add_argument("--project_templ", "-t", action="store_true",
                        default=False, help="enable projection of templates")
    args = parser.parse_args()

    if args.project:
        project(args.cfg_file)

    if args.project_templ:
        get_templates(args.cfg_file)

    if args.fit:
        fit(args.cfg_file)
