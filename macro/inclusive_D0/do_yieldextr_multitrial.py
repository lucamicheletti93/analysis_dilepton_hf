"""
Script to perform multi-trial for signal extraction as a function of pT
"""

import os
import sys
import argparse
import itertools
import numpy as np
import yaml
import pdg
from alive_progress import alive_bar
from fit_utils import perform_fit
import ROOT

def multi_trial(input_config):
    """
    """

    with open(input_config, "r") as f:
        cfg = yaml.safe_load(f)

    if cfg["hadron"] not in ["dzero"]:
        print(f"ERROR: {cfg['hadron']} not supported, exit")
        sys.exit()

    pdg_api = pdg.connect()

    pt_mins = cfg["pt_mins"]
    pt_maxs = cfg["pt_maxs"]
    infile_data = ROOT.TFile.Open(cfg["inputs"]["data"])
    infile_mc = ROOT.TFile.Open(cfg["inputs"]["mc"])
    for pt_min, pt_max in zip(pt_mins, pt_maxs):
        hist_data = infile_data.Get(f"hist_mass_pt{pt_min:.1f}_{pt_max:.1f}_nocutnp")
        hist_mc = infile_mc.Get(f"hist_signal_pt{pt_min:.1f}_{pt_max:.1f}")
        if not isinstance(hist_data, ROOT.TH1) or not isinstance(hist_mc, ROOT.TH1):
            print("ERROR: pT interval {} {} not found")
            sys.exit()
    infile_data.Close()
    infile_mc.Close()

    n_trials = len(cfg["sgn_funcs"]) * len(cfg["bkg_funcs"]) * len(cfg["mass_mins"]) * \
        len(cfg["mass_maxs"]) * len(cfg["rebin"]) * len(cfg["fix_frac_bkgcorr"]) * len(cfg["fix_tailpars_from_mc"])
    pt_limits = pt_mins.copy()
    pt_limits.append(pt_maxs[-1])
    pt_limits = np.array(pt_limits, np.float64)

    hist_rawyield, hist_mean, hist_sigma, hist_chi2 = [], [], [], []
    for itrial in range(n_trials):
        hist_rawyield.append(ROOT.TH1D(
            f"hist_rawyield_{itrial}", ";#it{p}_{T} (GeV/#it{c}); raw yield",
            len(pt_mins), pt_limits))
        hist_sigma.append(ROOT.TH1D(
            f"hist_sigma_{itrial}", ";#it{p}_{T} (GeV/#it{c}); #sigma (GeV/#it{c}^{2})",
            len(pt_mins), pt_limits))
        hist_mean.append(ROOT.TH1D(
            f"hist_mean_{itrial}", ";#it{p}_{T} (GeV/#it{c}); #mu (GeV/#it{c}^{2})",
            len(pt_mins), pt_limits))
        hist_chi2.append(ROOT.TH1D(
            f"hist_chi2_{itrial}", ";#it{p}_{T} (GeV/#it{c}); #chi^{2}/ndf",
            len(pt_mins), pt_limits))

    for ipt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
        with alive_bar(n_trials, title=f"PT BIN {ipt} [{pt_min:.1f}-{pt_max:.1f}]: Testing trial ") as bar:
            for itrial, pars in enumerate(itertools.product(cfg["sgn_funcs"],
                                                            cfg["bkg_funcs"],
                                                            cfg["mass_mins"],
                                                            cfg["mass_maxs"],
                                                            cfg["rebin"],
                                                            cfg["fix_frac_bkgcorr"],
                                                            cfg["fix_tailpars_from_mc"])):
                fitter, _ = perform_fit(cfg["inputs"]["data"],
                                        f"hist_mass_pt{pt_min:.1f}_{pt_max:.1f}_nocutnp",
                                        f"dzerofit_pt{pt_min:.1f}_{pt_max:.1f}_{itrial}",
                                        pars[0],
                                        pars[1],
                                        pars[2],
                                        pars[3],
                                        pdg_api.get_particle_by_mcid(421).mass,
                                        cfg["hadron"],
                                        cfg["inputs"]["mc"],
                                        f"hist_templ_refl_pt{pt_min:.1f}_{pt_max:.1f}",
                                        f"hist_signal_pt{pt_min:.1f}_{pt_max:.1f}",
                                        True,
                                        pars[6],
                                        {},
                                        pars[4],
                                        pars[5])
                if fitter.get_fit_result.converged:
                    rawyield = fitter.get_raw_yield(0)
                    sigma = fitter.get_signal_parameter(0, "sigma")
                    mean = fitter.get_mass(0)
                    hist_rawyield[itrial].SetBinContent(ipt+1, rawyield[0])
                    hist_rawyield[itrial].SetBinError(ipt+1, rawyield[1])
                    hist_sigma[itrial].SetBinContent(ipt+1, sigma[0])
                    hist_sigma[itrial].SetBinError(ipt+1, sigma[1])
                    hist_mean[itrial].SetBinContent(ipt+1, mean[0])
                    hist_mean[itrial].SetBinError(ipt+1, mean[1])
                    hist_chi2[itrial].SetBinContent(ipt+1, fitter.get_chi2_ndf())
                    hist_chi2[itrial].SetBinError(ipt+1, 1.e-20)
                bar()

    outfile = ROOT.TFile(
        os.path.join(cfg["output"]["directory"], cfg["output"]["file"]), "recreate")
    for itrial in range(n_trials):
        hist_rawyield[itrial].Write()
        hist_sigma[itrial].Write()
        hist_mean[itrial].Write()
        hist_chi2[itrial].Write()
    outfile.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--cfg_file", "-c", metavar="text",
                        default="config.yml", help="config file with multitrial config")
    args = parser.parse_args()

    multi_trial(args.cfg_file)
