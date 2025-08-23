"""
Script for the production of the preliminary plots for D-JPsi measurement
"""

import argparse
import numpy as np
import yaml
import ROOT


def set_style(hist, color, markerstyle=ROOT.kFullCircle, marksersize=0, fillstyle=0, linewidth=2, alpha=1):
    """
    Simple method to set the style
    """

    if isinstance(hist, ROOT.TH1):
        hist.SetDirectory(0)

    hist.SetLineWidth(linewidth)
    hist.SetMarkerStyle(markerstyle)
    hist.SetMarkerSize(marksersize)
    hist.SetLineColorAlpha(color, alpha)
    hist.SetMarkerColorAlpha(color, alpha)
    hist.SetFillStyle(fillstyle)

    if fillstyle > 0:
        hist.SetFillColorAlpha(color, alpha)

def sum_in_quadrature(array):
    """
    Helper method to compute sum in quadrature
    """
    
    sum_sq = 0
    for element in array:
        sum_sq += element**2

    return np.sqrt(sum_sq)

def integrate_hist(hist, unc_corr=False):
    """
    Method for fast integration of differential cross sections stored in a histogram
    """

    xsec, xsec_unc = 0., 0.
    for ibin in range(1, hist.GetNbinsX()+1):
        xsec += hist.GetBinContent(ibin) * hist.GetBinWidth(ibin)
        if not unc_corr:
            xsec_unc += hist.GetBinError(ibin)**2 * hist.GetBinWidth(ibin)**2
        else:
            xsec_unc += hist.GetBinError(ibin) * hist.GetBinWidth(ibin)
    if not unc_corr:
        xsec_unc = np.sqrt(xsec_unc)

    return xsec, xsec_unc


def integrate_graph(graph, unc_corr=True):
    """
    Method for fast integration of differential cross sections stored in a histogram
    """

    xsec, xsec_unc_low, xsec_unc_high = 0., 0., 0.
    for ibin in range(graph.GetN()):
        delta = (graph.GetErrorXlow(ibin) + graph.GetErrorXhigh(ibin))
        xsec += graph.GetPointY(ibin) * delta
        if not unc_corr:
            xsec_unc_low += graph.GetErrorYlow(ibin)**2 * delta**2
            xsec_unc_high += graph.GetErrorYhigh(ibin)**2 * delta**2
        else:
            xsec_unc_low += graph.GetErrorYlow(ibin) * delta
            xsec_unc_high += graph.GetErrorYhigh(ibin) * delta

    if not unc_corr:
        xsec_unc_low = np.sqrt(xsec_unc_low)
        xsec_unc_high = np.sqrt(xsec_unc_high)

    return xsec, xsec_unc_low, xsec_unc_high


def plot(sigma_eff_helac):
    """
    Main function for plots
    """

    ROOT.gStyle.SetTitleSize(0.035, "xy")
    ROOT.gStyle.SetLabelSize(0.035, "xy")
    ROOT.gStyle.SetPadRightMargin(0.035)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadTopMargin(0.035)
    ROOT.gStyle.SetTitleOffset(1.2, "x")
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

    # Dzero
    infile_dzero = ROOT.TFile.Open("../data_shared/dzero_xsec_pp13dot6TeV_y06.root")
    hist_dzero_xsec = infile_dzero.Get("hist_xsec_ptint_p")
    hist_dzero_xsec.SetDirectory(0)
    infile_dzero.Close()
    
    # JPsi
    infile_jpsi = ROOT.TFile.Open("../data_shared/jpsi_xsec_pp13dot6TeV.root")
    hist_jpsi_xsec = infile_jpsi.Get("histStatXsecInt")
    hist_jpsi_xsec.SetDirectory(0)
    hist_jpsi_xsec.Scale(1.5) # from dsigma/dy to sigma
    infile_jpsi.Close()

    # JPsi/D0
    infile_jpsiszero = ROOT.TFile.Open("../data_shared/jpsidzero_xsec_pp13dot6TeV.root")
    hist_jpsidzero_xsec = infile_jpsiszero.Get("histStatXsecInt")
    hist_jpsidzero_xsec.SetDirectory(0)
    hist_jpsidzero_ratio = infile_jpsiszero.Get("histStatRatioInt")
    hist_jpsidzero_ratio.SetDirectory(0)
    infile_jpsiszero.Close()

    # JPsi/D0
    # hist_jpsi_dzero_ratio = hist_jpsi_xsec.Clone("hist_jpsi_dzero_ratio")
    # hist_jpsi_dzero_ratio.SetDirectory(0)
    # hist_jpsi_dzero_ratio.Divide(hist_dzero_xsec)

    # load systematics
    with open("syst_db.yml", "r") as file_sys:
        cfg_sys = yaml.safe_load(file_sys)
    syst_dzero_tot = sum_in_quadrature([cfg_sys["dzero"]["yield_extraction"],
                                        cfg_sys["dzero"]["sel_eff"],
                                        cfg_sys["dzero"]["trk_eff"],
                                        cfg_sys["dzero"]["match_eff"],
                                        cfg_sys["dzero"]["frac"],
                                        cfg_sys["dzero"]["br"],
                                        cfg_sys["dzero"]["lumi"]])
    syst_jpsi_tot = sum_in_quadrature([cfg_sys["jpsi"]["yield_extraction"],
                                       cfg_sys["jpsi"]["trk_eff"],
                                       cfg_sys["jpsi"]["sel_eff"],
                                       cfg_sys["jpsi"]["match_eff"],
                                       cfg_sys["jpsi"]["br"],
                                       cfg_sys["jpsi"]["lumi"]])
    syst_jpsidmes_tot = sum_in_quadrature([cfg_sys["dzerojpsi"]["yield_extraction"],
                                           cfg_sys["dzero"]["sel_eff"],
                                           cfg_sys["dzero"]["trk_eff"],
                                           cfg_sys["dzero"]["match_eff"],
                                           cfg_sys["dzero"]["frac"],
                                           cfg_sys["dzero"]["br"],
                                           cfg_sys["dzero"]["lumi"],
                                           cfg_sys["jpsi"]["trk_eff"],
                                           cfg_sys["jpsi"]["sel_eff"],
                                           cfg_sys["jpsi"]["match_eff"],
                                           cfg_sys["jpsi"]["br"]])
    syst_jpsi_dzero_singlehad = sum_in_quadrature([cfg_sys["dzero"]["yield_extraction"],
                                                   cfg_sys["dzero"]["sel_eff"],
                                                   cfg_sys["dzero"]["trk_eff"],
                                                   cfg_sys["dzero"]["match_eff"],
                                                   cfg_sys["dzero"]["frac"],
                                                   cfg_sys["dzero"]["br"],
                                                   cfg_sys["jpsi"]["sel_eff"],
                                                   cfg_sys["jpsi"]["trk_eff"],
                                                   cfg_sys["jpsi"]["match_eff"],
                                                   cfg_sys["jpsi"]["br"]])
    syst_jpsi_dzero_ratio = sum_in_quadrature([cfg_sys["dzerojpsi"]["yield_extraction"],
                                               cfg_sys["dzero"]["lumi"]])

    #JPsi-D0
    graph_xsec_stat_all = ROOT.TGraphAsymmErrors(1)
    graph_xsec_syst_all = ROOT.TGraphAsymmErrors(1)
    graph_xsec_stat_all.SetName("graph_xsec_stat_all")
    graph_xsec_syst_all.SetName("graph_xsec_stat_all")
    set_style(graph_xsec_stat_all, ROOT.kBlack, ROOT.kFullCircle, 1.5, 0, 2, 1)
    set_style(graph_xsec_syst_all, ROOT.kBlack, ROOT.kFullCircle, 1.5, 0, 2, 1)
    xsec_dzero = hist_dzero_xsec.GetBinContent(1)
    xsec_jpsi = hist_jpsi_xsec.GetBinContent(1)
    graph_xsec_stat_all.SetPoint(0, xsec_dzero, 4.) # D meson in fiducial acceptance
    graph_xsec_stat_all.SetPoint(1, xsec_jpsi, 3.) # J/Psi meson in fiducial acceptance
    graph_xsec_syst_all.SetPoint(0, xsec_dzero, 4.) # D meson in fiducial acceptance
    graph_xsec_syst_all.SetPoint(1, xsec_jpsi, 3.) # J/Psi meson in fiducial acceptance
    #xsec_jpsidmes = 0.282189 # TODO: replace me with a value from file
    #xsec_unc_jpsidmes = xsec_jpsidmes * 0.07 # TODO: replace me with a value from file
    xsec_jpsidmes = hist_jpsidzero_xsec.GetBinContent(1)
    xsec_unc_jpsidmes = hist_jpsidzero_xsec.GetBinError(1)
    graph_xsec_stat_all.SetPoint(2, xsec_jpsidmes, 2.)
    graph_xsec_syst_all.SetPoint(2, xsec_jpsidmes, 2.)
    graph_xsec_stat_all.SetPointError(0, hist_dzero_xsec.GetBinError(1), hist_dzero_xsec.GetBinError(1), 0., 0.) # D meson in fiducial acceptance
    graph_xsec_stat_all.SetPointError(1, hist_jpsi_xsec.GetBinError(1), hist_jpsi_xsec.GetBinError(1), 0., 0.) # J/Psi meson in fiducial acceptance
    graph_xsec_stat_all.SetPointError(2, xsec_unc_jpsidmes, xsec_unc_jpsidmes, 0., 0.) # TODO: replace me with a value from file
    graph_xsec_syst_all.SetPointError(0, syst_dzero_tot * xsec_dzero, syst_dzero_tot * xsec_dzero, 0.1, 0.1) # D meson in fiducial acceptance
    graph_xsec_syst_all.SetPointError(1, syst_jpsi_tot * xsec_jpsi, syst_jpsi_tot * xsec_jpsi, 0.1, 0.1) # D meson in fiducial acceptance
    graph_xsec_syst_all.SetPointError(2, syst_jpsidmes_tot * xsec_jpsidmes, syst_jpsidmes_tot * xsec_jpsidmes, 0.1, 0.1) # TODO: replace me with a value from file

    graph_xsec_deltay_stat = ROOT.TGraphAsymmErrors(1)
    graph_xsec_deltay_stat.SetName("graph_xsec_deltay_stat")
    graph_xsec_deltay_syst = ROOT.TGraphAsymmErrors(1)
    graph_xsec_deltay_syst.SetName("graph_xsec_deltay_syst")
    set_style(graph_xsec_deltay_stat, ROOT.kBlack, ROOT.kFullCircle, 1.5, 0, 2, 1)
    set_style(graph_xsec_deltay_syst, ROOT.kBlack, ROOT.kFullCircle, 1.5, 0, 2, 1)
    graph_xsec_deltay_stat.SetPoint(0, (4.6 + 1.9) / 2, xsec_jpsidmes / (4.6 - 1.9))
    graph_xsec_deltay_syst.SetPoint(0, (4.6 + 1.9) / 2, xsec_jpsidmes / (4.6 - 1.9))
    graph_xsec_deltay_stat.SetPointError(0, (4.6 - 1.9) / 2, (4.6 - 1.9) / 2,
                                         xsec_unc_jpsidmes / (4.6 - 1.9), xsec_unc_jpsidmes / (4.6 - 1.9))
    graph_xsec_deltay_syst.SetPointError(0, (4.6 - 1.9) / 4, (4.6 - 1.9) / 4,
                                         syst_jpsidmes_tot * xsec_jpsidmes / (4.6 - 1.9), syst_jpsidmes_tot * xsec_jpsidmes / (4.6 - 1.9))

    graph_ratio_singlehad_stat = ROOT.TGraphAsymmErrors(1)
    graph_ratio_singlehad_stat.SetName("graph_ratio_singlehad_stat")
    graph_ratio_singlehad_syst = ROOT.TGraphAsymmErrors(1)
    graph_ratio_singlehad_syst.SetName("graph_ratio_singlehad_syst")
    set_style(graph_ratio_singlehad_stat, ROOT.kBlack, ROOT.kFullCircle, 1.5, 0, 2, 1)
    set_style(graph_ratio_singlehad_syst, ROOT.kBlack, ROOT.kFullCircle, 1.5, 0, 2, 1)
    ratio = hist_jpsi_xsec.GetBinContent(1)/hist_dzero_xsec.GetBinContent(1)
    graph_ratio_singlehad_stat.SetPoint(0, ratio, 2.)
    graph_ratio_singlehad_syst.SetPoint(0, ratio, 2.)
    stat_unc = np.sqrt((
        hist_jpsi_xsec.GetBinError(1)/hist_jpsi_xsec.GetBinContent(1))**2 + (hist_dzero_xsec.GetBinError(1)/hist_dzero_xsec.GetBinContent(1))**2) * ratio
    graph_ratio_singlehad_stat.SetPointError(0, stat_unc, stat_unc, 0., 0.)
    graph_ratio_singlehad_syst.SetPointError(0, syst_jpsi_dzero_singlehad * ratio,
                                                syst_jpsi_dzero_singlehad * ratio, 0.1, 0.1)

    graph_ratio_stat = ROOT.TGraphAsymmErrors(1)
    graph_ratio_stat.SetName("graph_ratio_stat")
    graph_ratio_syst = ROOT.TGraphAsymmErrors(1)
    graph_ratio_syst.SetName("graph_ratio_syst")
    set_style(graph_ratio_stat, ROOT.kBlack, ROOT.kFullCircle, 1.5, 0, 2, 1)
    set_style(graph_ratio_syst, ROOT.kBlack, ROOT.kFullCircle, 1.5, 0, 2, 1)
    #ratio = 37.481 # TODO: replace me with a value from file
    ratio = hist_jpsidzero_ratio.GetBinContent(1) / 1000.
    graph_ratio_stat.SetPoint(0, ratio, 2.)
    graph_ratio_syst.SetPoint(0, ratio, 2.)
    #stat_unc = 0.07 * ratio # TODO: replace me with a value from file
    stat_unc = hist_jpsidzero_ratio.GetBinError(1) / 1000.
    graph_ratio_stat.SetPointError(0, stat_unc, stat_unc, 0., 0.)
    graph_ratio_syst.SetPointError(0, syst_jpsi_dzero_ratio * ratio, syst_jpsi_dzero_ratio * ratio, 0.1, 0.1)

    # LHCb
    ratio_lhcb = 14.9 # mb
    syst_low_ratio_lhcb = np.sqrt(1.1**2 + 3.1**2)
    syst_high_ratio_lhcb = np.sqrt(1.1**2 + 2.3**2)
    stat_ratio_lhcb = 0.4
    graph_ratio_stat.SetPoint(1, ratio_lhcb, 1.)
    graph_ratio_syst.SetPoint(1, ratio_lhcb, 1.)
    graph_ratio_stat.SetPointError(1, stat_ratio_lhcb, stat_ratio_lhcb, 0., 0.)
    graph_ratio_syst.SetPointError(1, syst_low_ratio_lhcb, syst_high_ratio_lhcb, 0.1, 0.1)

    # Simulations
    infile_sim = ROOT.TFile.Open("../sim/predictions_pythia8_jpsidmes_kMonash_kSoftQCD.root")
    hist_pythia_dmes_xsec = infile_sim.Get("hist_dmesprompt_ptint")
    hist_pythia_jpsi_fwd_xsec = infile_sim.Get("hist_jpsiincl_ptint")
    hist_pythia_dmes_jpsi_fwd_xsec_tot = infile_sim.Get("fwd_cuts/hist_deltay_tot_jpsiincl_dmesprompt")
    hist_pythia_dmes_jpsi_fwd_xsec_dps = infile_sim.Get("fwd_cuts/hist_deltay_dps_jpsiincl_dmesprompt")
    hist_pythia_dmes_jpsi_fwd_xsec_sps = infile_sim.Get("fwd_cuts/hist_deltay_sps_jpsiincl_dmesprompt")
    hist_pythia_dmes_jpsi_mid_xsec_tot = infile_sim.Get("mid_cuts/hist_deltay_tot_jpsiincl_dmesprompt")
    hist_pythia_dmes_jpsi_mid_xsec_dps = infile_sim.Get("mid_cuts/hist_deltay_dps_jpsiincl_dmesprompt")
    hist_pythia_dmes_jpsi_mid_xsec_sps = infile_sim.Get("mid_cuts/hist_deltay_sps_jpsiincl_dmesprompt")
    graph_helac_dmes_jpsi_fwd_xsec_dps = infile_sim.Get("fwd_cuts/graph_helac_deltay_dps_jpsiincl_dmesincl")
    graph_helac_dmes_jpsi_fwd_xsec_sps = infile_sim.Get("fwd_cuts/graph_helac_deltay_sps_jpsiprompt_dmesprompt")
    graph_helac_dmes_jpsi_mid_xsec_dps = infile_sim.Get("mid_cuts/graph_helac_deltay_dps_jpsiincl_dmesincl")
    graph_helac_dmes_jpsi_mid_xsec_sps = infile_sim.Get("mid_cuts/graph_helac_deltay_sps_jpsiprompt_dmesprompt")
    hist_pythia_dmes_jpsi_ratio_tot = infile_sim.Get("hist_ratio_tot_jpsiincl_dmesprompt")
    hist_pythia_dmes_jpsi_ratio_dps = infile_sim.Get("hist_ratio_dps_jpsiincl_dmesprompt")
    hist_pythia_dmes_jpsi_ratio_sps = infile_sim.Get("hist_ratio_sps_jpsiincl_dmesprompt")
    hist_pythia_dmes_jpsi_ratiosinglehad_tot = infile_sim.Get("hist_ratio_singlehad_jpsiincl_dmesprompt")
    hist_pythia_dmes_xsec.SetDirectory(0)
    hist_pythia_jpsi_fwd_xsec.SetDirectory(0)
    hist_pythia_dmes_jpsi_fwd_xsec_tot.SetDirectory(0)
    hist_pythia_dmes_jpsi_fwd_xsec_dps.SetDirectory(0)
    hist_pythia_dmes_jpsi_fwd_xsec_sps.SetDirectory(0)
    hist_pythia_dmes_jpsi_mid_xsec_tot.SetDirectory(0)
    hist_pythia_dmes_jpsi_mid_xsec_dps.SetDirectory(0)
    hist_pythia_dmes_jpsi_mid_xsec_sps.SetDirectory(0)
    hist_pythia_dmes_jpsi_ratio_tot.SetDirectory(0)
    hist_pythia_dmes_jpsi_ratio_dps.SetDirectory(0)
    hist_pythia_dmes_jpsi_ratio_sps.SetDirectory(0)
    hist_pythia_dmes_jpsi_ratiosinglehad_tot.SetDirectory(0)
    infile_sim.Close()

    infile_sim_7tev = ROOT.TFile.Open("../sim/predictions_pythia8_jpsidmes_kMonash_kSoftQCD_7Tev.root")
    hist_pythia_dmes_jpsi_ratio_tot_7tev = infile_sim_7tev.Get("hist_ratio_tot_jpsiincl_dmesprompt")
    hist_pythia_dmes_jpsi_ratio_dps_7tev = infile_sim_7tev.Get("hist_ratio_dps_jpsiincl_dmesprompt")
    hist_pythia_dmes_jpsi_ratio_sps_7tev = infile_sim_7tev.Get("hist_ratio_sps_jpsiincl_dmesprompt")
    hist_pythia_dmes_jpsi_ratiosinglehad_tot_7tev = infile_sim_7tev.Get("hist_ratio_singlehad_jpsiprompt_dmesprompt")
    hist_pythia_dmes_jpsi_ratio_tot_7tev.SetDirectory(0)
    hist_pythia_dmes_jpsi_ratio_dps_7tev.SetDirectory(0)
    hist_pythia_dmes_jpsi_ratio_sps_7tev.SetDirectory(0)
    hist_pythia_dmes_jpsi_ratiosinglehad_tot_7tev.SetDirectory(0)
    infile_sim_7tev.Close()

    # PYTHIA
    graph_pythia_xsec_all = ROOT.TGraphAsymmErrors(1)
    graph_pythia_xsec_all.SetName("graph_pythia_xsec_all")
    set_style(graph_pythia_xsec_all, ROOT.kGray+1, 1, 0, 0, 4, 0.6)
    graph_pythia_xsec_all.SetPoint(0, hist_pythia_dmes_xsec.GetBinContent(2), 4.) # D meson in fiducial acceptance
    graph_pythia_xsec_all.SetPoint(1, hist_pythia_jpsi_fwd_xsec.GetBinContent(2), 3.) # J/Psi meson in fiducial acceptance
    graph_pythia_xsec_all.SetPointError(0, hist_pythia_dmes_xsec.GetBinError(2),
                                        hist_pythia_dmes_xsec.GetBinError(2), 0.25, 0.25) # D meson in fiducial acceptance
    graph_pythia_xsec_all.SetPointError(1, hist_pythia_jpsi_fwd_xsec.GetBinError(2),
                                        hist_pythia_jpsi_fwd_xsec.GetBinError(2), 0.25, 0.25) # J/Psi meson in fiducial acceptance
    xsec, xsec_unc = integrate_hist(hist_pythia_dmes_jpsi_fwd_xsec_tot)
    graph_pythia_xsec_all.SetPoint(2, xsec, 2.) # D meson in fiducial acceptance
    graph_pythia_xsec_all.SetPointError(2, xsec_unc, xsec_unc, 0.25, 0.25) # D meson in fiducial acceptance

    graph_pythia_xsec_dps = ROOT.TGraphAsymmErrors(1)
    graph_pythia_xsec_dps.SetName("graph_pythia_xsec_dps")
    set_style(graph_pythia_xsec_dps, ROOT.kBlue+2, 1, 0, 0, 4, 0.6)
    xsec, xsec_unc = integrate_hist(hist_pythia_dmes_jpsi_fwd_xsec_dps)
    graph_pythia_xsec_dps.SetPoint(0, xsec, 2.) # D meson in fiducial acceptance
    graph_pythia_xsec_dps.SetPointError(0, xsec_unc, xsec_unc, 0.25, 0.25) # D meson in fiducial acceptance

    graph_pythia_xsec_sps = ROOT.TGraphAsymmErrors(1)
    graph_pythia_xsec_sps.SetName("graph_pythia_xsec_dps")
    set_style(graph_pythia_xsec_sps, ROOT.kRed+1, 1, 0, 0, 4, 0.6)
    xsec, xsec_unc = integrate_hist(hist_pythia_dmes_jpsi_fwd_xsec_sps)
    graph_pythia_xsec_sps.SetPoint(0, xsec, 2.) # D meson in fiducial acceptance
    graph_pythia_xsec_sps.SetPointError(0, xsec_unc, xsec_unc, 0.25, 0.25) # D meson in fiducial acceptance

    graph_pythia_xsec_deltay_tot = ROOT.TGraphAsymmErrors(1)
    graph_pythia_xsec_deltay_tot.SetName("graph_pythia_xsec_deltay_tot")
    set_style(graph_pythia_xsec_deltay_tot, ROOT.kGray+1, 1, 0, 0, 4, 0.6)
    xsec, xsec_unc = integrate_hist(hist_pythia_dmes_jpsi_mid_xsec_tot)
    deltay = 1.5+1.5
    graph_pythia_xsec_deltay_tot.SetPoint(0, 0., xsec / deltay)
    graph_pythia_xsec_deltay_tot.SetPointError(0, 1.5, 1.5, xsec_unc / deltay, xsec_unc / deltay)
    xsec, xsec_unc = integrate_hist(hist_pythia_dmes_jpsi_fwd_xsec_tot)
    deltay = 4.6-1.9
    graph_pythia_xsec_deltay_tot.SetPoint(1, 3.25, xsec / deltay)
    graph_pythia_xsec_deltay_tot.SetPointError(1, 1.35, 1.35, xsec_unc / deltay, xsec_unc / deltay)

    graph_pythia_xsec_deltay_dps = ROOT.TGraphAsymmErrors(1)
    graph_pythia_xsec_deltay_dps.SetName("graph_pythia_xsec_deltay_dps")
    set_style(graph_pythia_xsec_deltay_dps, ROOT.kBlue+2, 1, 0, 0, 4, 0.6)
    xsec, xsec_unc = integrate_hist(hist_pythia_dmes_jpsi_mid_xsec_dps)
    deltay = 1.5+1.5
    graph_pythia_xsec_deltay_dps.SetPoint(0, 0., xsec / deltay)
    graph_pythia_xsec_deltay_dps.SetPointError(0, 1.5, 1.5, xsec_unc / deltay, xsec_unc / deltay)
    xsec, xsec_unc = integrate_hist(hist_pythia_dmes_jpsi_fwd_xsec_dps)
    deltay = 4.6-1.9
    graph_pythia_xsec_deltay_dps.SetPoint(1, 3.25, xsec / deltay)
    graph_pythia_xsec_deltay_dps.SetPointError(1, 1.35, 1.35, xsec_unc / deltay, xsec_unc / deltay)

    graph_pythia_xsec_deltay_sps = ROOT.TGraphAsymmErrors(1)
    graph_pythia_xsec_deltay_sps.SetName("graph_pythia_xsec_deltay_sps")
    set_style(graph_pythia_xsec_deltay_sps, ROOT.kRed+1, 1, 0, 0, 4, 0.6)
    xsec, xsec_unc = integrate_hist(hist_pythia_dmes_jpsi_mid_xsec_sps)
    deltay = 1.5+1.5
    graph_pythia_xsec_deltay_sps.SetPoint(0, 0., xsec / deltay)
    graph_pythia_xsec_deltay_sps.SetPointError(0, 1.5, 1.5, xsec_unc / deltay, xsec_unc / deltay)
    xsec, xsec_unc = integrate_hist(hist_pythia_dmes_jpsi_fwd_xsec_sps)
    deltay = 4.6-1.9
    graph_pythia_xsec_deltay_sps.SetPoint(1, 3.25, xsec / deltay)
    graph_pythia_xsec_deltay_sps.SetPointError(1, 1.35, 1.35, xsec_unc / deltay, xsec_unc / deltay)

    graph_pythia_ratio_tot = ROOT.TGraphAsymmErrors(1)
    graph_pythia_ratio_tot.SetName("graph_pythia_ratio_tot")
    set_style(graph_pythia_ratio_tot, ROOT.kGray+1, 1, 0, 1000, 4, 0.6)
    graph_pythia_ratio_tot.SetPoint(0, hist_pythia_dmes_jpsi_ratio_tot.GetBinContent(2)/1000, 2.) #mb
    graph_pythia_ratio_tot.SetPointError(0, hist_pythia_dmes_jpsi_ratio_tot.GetBinError(2)/1000,
                                         hist_pythia_dmes_jpsi_ratio_tot.GetBinError(2)/1000, 0.25, 0.25)
    graph_pythia_ratio_tot.SetPoint(1, hist_pythia_dmes_jpsi_ratio_tot_7tev.GetBinContent(2)/1000, 1.) #mb
    graph_pythia_ratio_tot.SetPointError(1, hist_pythia_dmes_jpsi_ratio_tot_7tev.GetBinError(2)/1000,
                                         hist_pythia_dmes_jpsi_ratio_tot_7tev.GetBinError(2)/1000, 0.25, 0.25)

    graph_pythia_ratio_sps = ROOT.TGraphAsymmErrors(1)
    graph_pythia_ratio_sps.SetName("graph_pythia_ratio_sps")
    set_style(graph_pythia_ratio_sps, ROOT.kRed+1, 1, 0, 1000, 4, 0.6)
    graph_pythia_ratio_sps.SetPoint(0, hist_pythia_dmes_jpsi_ratio_sps.GetBinContent(2)/1000, 2.) #mb
    graph_pythia_ratio_sps.SetPointError(0, hist_pythia_dmes_jpsi_ratio_sps.GetBinError(2)/1000,
                                         hist_pythia_dmes_jpsi_ratio_sps.GetBinError(2)/1000, 0.25, 0.25)
    graph_pythia_ratio_sps.SetPoint(1, hist_pythia_dmes_jpsi_ratio_sps_7tev.GetBinContent(2)/1000, 1.) #mb
    graph_pythia_ratio_sps.SetPointError(1, hist_pythia_dmes_jpsi_ratio_sps_7tev.GetBinError(2)/1000,
                                         hist_pythia_dmes_jpsi_ratio_sps_7tev.GetBinError(2)/1000, 0.25, 0.25)

    graph_pythia_ratio_dps = ROOT.TGraphAsymmErrors(1)
    graph_pythia_ratio_dps.SetName("graph_pythia_ratio_dps")
    set_style(graph_pythia_ratio_dps, ROOT.kBlue+1, 1, 0, 1000, 4, 0.6)
    graph_pythia_ratio_dps.SetPoint(0, hist_pythia_dmes_jpsi_ratio_dps.GetBinContent(2)/1000, 2.) #mb
    graph_pythia_ratio_dps.SetPointError(0, hist_pythia_dmes_jpsi_ratio_dps.GetBinError(2)/1000,
                                         hist_pythia_dmes_jpsi_ratio_dps.GetBinError(2)/1000, 0.25, 0.25)
    graph_pythia_ratio_dps.SetPoint(1, hist_pythia_dmes_jpsi_ratio_dps_7tev.GetBinContent(2)/1000, 1.) #mb
    graph_pythia_ratio_dps.SetPointError(1, hist_pythia_dmes_jpsi_ratio_dps_7tev.GetBinError(2)/1000,
                                         hist_pythia_dmes_jpsi_ratio_dps_7tev.GetBinError(2)/1000, 0.25, 0.25)

    graph_pythia_ratio_singlehad = ROOT.TGraphAsymmErrors(1)
    graph_pythia_ratio_singlehad.SetName("graph_pythia_ratio_singlehad")
    set_style(graph_pythia_ratio_singlehad, ROOT.kGray+1, 1, 0, 1000, 0, 0.6)
    graph_pythia_ratio_singlehad.SetPoint(0, hist_pythia_dmes_jpsi_ratiosinglehad_tot.GetBinContent(2), 2.)
    graph_pythia_ratio_singlehad.SetPointError(0, hist_pythia_dmes_jpsi_ratiosinglehad_tot.GetBinError(2),
                                               hist_pythia_dmes_jpsi_ratiosinglehad_tot.GetBinError(2), 0.25, 0.25)
    graph_pythia_ratio_singlehad.SetPoint(1, hist_pythia_dmes_jpsi_ratiosinglehad_tot_7tev.GetBinContent(2), 1.)
    graph_pythia_ratio_singlehad.SetPointError(1, hist_pythia_dmes_jpsi_ratiosinglehad_tot_7tev.GetBinError(2),
                                               hist_pythia_dmes_jpsi_ratiosinglehad_tot_7tev.GetBinError(2), 0.25, 0.25)

    # Helac ONIA
    xsec_mid_dps, xsec_unc_mid_dps_low, xsec_unc_mid_dps_high = integrate_graph(graph_helac_dmes_jpsi_mid_xsec_dps)
    xsec_mid_dps = xsec_mid_dps / sigma_eff_helac * 15
    xsec_unc_mid_dps_low = xsec_unc_mid_dps_low / sigma_eff_helac * 15
    xsec_unc_mid_dps_high = xsec_unc_mid_dps_high / sigma_eff_helac * 15
    xsec_fwd_dps, xsec_unc_fwd_dps_low, xsec_unc_fwd_dps_high = integrate_graph(graph_helac_dmes_jpsi_fwd_xsec_dps)
    xsec_fwd_dps = xsec_fwd_dps / sigma_eff_helac * 15
    xsec_unc_fwd_dps_low = xsec_unc_fwd_dps_low / sigma_eff_helac * 15
    xsec_unc_fwd_dps_high = xsec_unc_fwd_dps_high / sigma_eff_helac * 15
    xsec_mid_sps, xsec_unc_mid_sps_low, xsec_unc_mid_sps_high = integrate_graph(graph_helac_dmes_jpsi_mid_xsec_sps)
    xsec_mid_tot = xsec_mid_dps + xsec_mid_sps
    xsec_unc_mid_tot_low = np.sqrt(xsec_unc_mid_dps_low**2 + xsec_unc_mid_sps_low**2)
    xsec_unc_mid_tot_high = np.sqrt(xsec_unc_mid_dps_high**2 + xsec_unc_mid_sps_high**2)
    xsec_fwd_sps, xsec_unc_fwd_sps_low, xsec_unc_fwd_sps_high = integrate_graph(graph_helac_dmes_jpsi_fwd_xsec_sps)
    xsec_fwd_tot = xsec_fwd_dps + xsec_fwd_sps
    xsec_unc_fwd_tot_low = np.sqrt(xsec_unc_fwd_dps_low**2 + xsec_unc_fwd_sps_low**2)
    xsec_unc_fwd_tot_high = np.sqrt(xsec_unc_fwd_dps_high**2 + xsec_unc_fwd_sps_high**2)

    graph_helac_xsec_dps = ROOT.TGraphAsymmErrors(1)
    graph_helac_xsec_dps.SetName("graph_helac_xsec_dps")
    set_style(graph_helac_xsec_dps, ROOT.kAzure+1, 1, 0, 1000, 4, 0.6)
    graph_helac_xsec_dps.SetPoint(0, xsec_fwd_dps, 2.) # D meson in fiducial acceptance
    graph_helac_xsec_dps.SetPointError(0, xsec_unc_fwd_dps_low, xsec_unc_fwd_dps_high, 0.25, 0.25) # D meson in fiducial acceptance

    graph_helac_xsec_sps = ROOT.TGraphAsymmErrors(1)
    graph_helac_xsec_sps.SetName("graph_helac_xsec_sps")
    set_style(graph_helac_xsec_sps, ROOT.kOrange+6, 1, 0, 3154, 2, 0.4)
    graph_helac_xsec_sps.SetPoint(0, xsec_fwd_sps, 2.) # D meson in fiducial acceptance
    graph_helac_xsec_sps.SetPointError(0, xsec_unc_fwd_sps_low, xsec_unc_fwd_sps_high, 0.25, 0.25) # D meson in fiducial acceptance

    graph_helac_xsec_all = ROOT.TGraphAsymmErrors(1)
    graph_helac_xsec_all.SetName("graph_helac_xsec_all")
    set_style(graph_helac_xsec_all, ROOT.kTeal-7, 1, 0, 3145, 2, 0.4)
    graph_helac_xsec_all.SetPoint(1, xsec_fwd_tot, 2.) # D meson in fiducial acceptance
    graph_helac_xsec_all.SetPointError(1, xsec_unc_fwd_tot_low, xsec_unc_fwd_tot_high, 0.25, 0.25) # D meson in fiducial acceptance

    graph_helac_xsec_deltay_dps = ROOT.TGraphAsymmErrors(1)
    graph_helac_xsec_deltay_dps.SetName("graph_helac_xsec_deltay_dps")
    set_style(graph_helac_xsec_deltay_dps, ROOT.kAzure+1, 1, 0, 1000, 4, 0.6)
    deltay = 1.5+1.5
    graph_helac_xsec_deltay_dps.SetPoint(0, 0., xsec_mid_dps / deltay)
    graph_helac_xsec_deltay_dps.SetPointError(0, 1.5, 1.5, xsec_unc_mid_dps_low / deltay, xsec_unc_mid_dps_high / deltay)
    deltay = 4.6-1.9
    graph_helac_xsec_deltay_dps.SetPoint(1, 3.25, xsec_fwd_dps / deltay)
    graph_helac_xsec_deltay_dps.SetPointError(1, 1.35, 1.35, xsec_unc_fwd_dps_low / deltay, xsec_unc_fwd_dps_high / deltay)

    graph_helac_xsec_deltay_sps = ROOT.TGraphAsymmErrors(1)
    graph_helac_xsec_deltay_sps.SetName("graph_helac_xsec_deltay_sps")
    set_style(graph_helac_xsec_deltay_sps, ROOT.kOrange+6, 1, 0, 3154, 2, 0.4)
    deltay = 1.5+1.5
    graph_helac_xsec_deltay_sps.SetPoint(0, 0., xsec_mid_sps / deltay)
    graph_helac_xsec_deltay_sps.SetPointError(0, 1.5, 1.5, xsec_unc_mid_sps_low / deltay, xsec_unc_mid_sps_high / deltay)
    deltay = 4.6-1.9
    graph_helac_xsec_deltay_sps.SetPoint(1, 3.25, xsec_fwd_sps / deltay)
    graph_helac_xsec_deltay_sps.SetPointError(1, 1.35, 1.35, xsec_unc_fwd_sps_low / deltay, xsec_unc_fwd_sps_high / deltay)

    graph_helac_xsec_deltay_tot = ROOT.TGraphAsymmErrors(1)
    graph_helac_xsec_deltay_tot.SetName("graph_helac_xsec_deltay_tot")
    set_style(graph_helac_xsec_deltay_tot, ROOT.kTeal-7, 1, 0, 3145, 2, 0.4)
    deltay = 1.5+1.5
    graph_helac_xsec_deltay_tot.SetPoint(0, 0., xsec_mid_tot / deltay)
    graph_helac_xsec_deltay_tot.SetPointError(0, 1.5, 1.5, xsec_unc_mid_tot_low / deltay, xsec_unc_mid_tot_high / deltay)
    deltay = 4.6-1.9
    graph_helac_xsec_deltay_tot.SetPoint(1, 3.25, xsec_fwd_tot / deltay)
    graph_helac_xsec_deltay_tot.SetPointError(1, 1.35, 1.35, xsec_unc_fwd_tot_low / deltay, xsec_unc_fwd_tot_high / deltay)

    # Plots
    lat_label = ROOT.TLatex()
    lat_label.SetNDC()
    lat_label.SetTextSize(0.045)
    lat_label.SetTextFont(42)
    lat_label.SetTextColor(ROOT.kBlack)
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextSize(0.035)
    lat.SetTextFont(42)
    lat.SetTextColor(ROOT.kBlack)

    canv_xsec = ROOT.TCanvas("canv_xsec", "", 1200, 1200)
    leg_pythia = ROOT.TLegend(0.32, 0.14, 0.6, 0.28)
    leg_pythia.SetTextSize(0.03)
    leg_pythia.SetBorderSize(0)
    leg_pythia.SetFillStyle(0)
    leg_pythia.SetHeader("#splitline{PYTHIA 8 Monash}{CharmoniumShower:all}")
    leg_pythia.SetNColumns(2)
    leg_pythia.AddEntry(graph_pythia_xsec_sps, "SPS", "l")
    leg_pythia.AddEntry(graph_pythia_xsec_dps, "DPS", "l")
    leg_pythia.AddEntry(graph_pythia_xsec_all, "Total", "l")
    leg_helac = ROOT.TLegend(0.65, 0.14, 0.95, 0.28)
    leg_helac.SetTextSize(0.03)
    leg_helac.SetBorderSize(0)
    leg_helac.SetFillStyle(0)
    leg_helac.SetHeader("HELAC-Onia")
    leg_helac.SetMargin(0.15)
    leg_helac.AddEntry(graph_helac_xsec_sps, "SPS LO + PYTHIA 8", "f")
    leg_helac.AddEntry(graph_helac_xsec_dps,
                       f"DPS #sigma_{{eff}} = {sigma_eff_helac:.0f} mb", "f")
    leg_helac.AddEntry(graph_helac_xsec_all, "Total", "f")

    canv_xsec.SetLeftMargin(0.3)
    frame_xsec = ROOT.TH2D("frame_xsec", ";#sigma (#mub);", 1000, 2.e-2, 4000., 5, 0.5, 5.5)
    frame_xsec.GetYaxis().SetBinLabel(
        4, "#splitline{Prompt D^{0}, |#it{y}| < 0.6}{#it{p}_{T} > 0.5 GeV/#it{c}}")
    frame_xsec.GetYaxis().SetBinLabel(
        3, "#splitline{Inclusive J/#Psi}{2.5 < y < 4, #it{p}_{T} > 0}")
    frame_xsec.GetYaxis().SetBinLabel(2, "#splitline{D^{0}#minusJ/#Psi}{1.9 < #Delta#it{y} < 4.6}")
    frame_xsec.GetYaxis().SetLabelSize(0.05)
    canv_xsec.SetLogx()
    frame_xsec.DrawCopy()
    lat_label.DrawLatex(0.33, 0.9, "ALICE Preliminary")
    lat.DrawLatex(0.33, 0.86, "pp collisions, #sqrt{s} = 13.6 TeV")
    lat.DrawLatex(0.33, 0.82, "#it{L}_{int} = 39.7 pb^{#minus1}")
    graph_pythia_xsec_all.Draw("pz")
    graph_pythia_xsec_dps.Draw("pz")
    graph_pythia_xsec_sps.Draw("pz")
    graph_helac_xsec_all.Draw("2")
    graph_helac_xsec_all_cent = graph_helac_xsec_all.Clone("graph_helac_xsec_all_cent")
    graph_helac_xsec_all_cent.SetLineWidth(4)
    graph_helac_xsec_all_cent.SetMarkerSize(0)
    for ibin in range(graph_helac_xsec_all_cent.GetN()):
        graph_helac_xsec_all_cent.SetPointError(ibin, 0., 0., graph_helac_xsec_all.GetErrorYlow(ibin),
                                                graph_helac_xsec_all.GetErrorYhigh(ibin))
    graph_helac_xsec_all_cent.Draw("pz")

    graph_helac_xsec_dps.Draw("pz")
    graph_helac_xsec_sps.Draw("2")
    graph_xsec_syst_all.Draw("2")
    graph_xsec_stat_all.Draw("pz")
    leg_pythia.Draw()
    leg_helac.Draw()
    canv_xsec.SaveAs(f"dzero_jpsi_xsection_pp13dot6TeV_HELACONIAsigmaeff{sigma_eff_helac:.0f}.pdf")

    canv_xsec_deltay = ROOT.TCanvas("canv_xsec_deltay", "", 1200, 1200)
    frame = canv_xsec_deltay.DrawFrame(-2., 0., 5., 1.15, ";#Delta#it{y};d#sigma/d#Delta#it{y} (#mub)")
    frame.GetYaxis().SetDecimals()
    lat_label.DrawLatex(0.18, 0.9, "ALICE Preliminary")
    lat.DrawLatex(0.18, 0.86, "pp collisions, #sqrt{s} = 13.6 TeV, #it{L}_{int} = 39.7 pb^{#minus1}")
    lat.DrawLatex(0.18, 0.81, "Prompt D^{0}, |#it{y}| < 0.6, #it{p}_{T} > 0.5 GeV/#it{c}")
    lat.DrawLatex(0.18, 0.77, "Inclusive J/#Psi, 2.5 < #it{y} < 4.0, #it{p}_{T} > 0")
    graph_pythia_xsec_deltay_tot.Draw("pz")
    graph_pythia_xsec_deltay_dps.Draw("pz")
    graph_pythia_xsec_deltay_sps.Draw("pz")
    graph_helac_xsec_deltay_sps.Draw("2")
    graph_helac_xsec_deltay_sps_cent = graph_helac_xsec_deltay_sps.Clone("graph_helac_xsec_deltay_sps")
    graph_helac_xsec_deltay_sps_cent.SetLineWidth(4)
    graph_helac_xsec_deltay_sps_cent.SetMarkerSize(0)
    for ibin in range(graph_helac_xsec_deltay_sps_cent.GetN()):
        graph_helac_xsec_deltay_sps_cent.SetPointError(ibin, graph_helac_xsec_deltay_sps_cent.GetErrorXlow(ibin),
                                                       graph_helac_xsec_deltay_sps_cent.GetErrorXhigh(ibin), 0., 0.)
    graph_helac_xsec_deltay_sps_cent.Draw("pz")
    graph_helac_xsec_deltay_dps.Draw("pz")
    graph_helac_xsec_deltay_tot.Draw("2")
    graph_helac_xsec_deltay_tot_cent = graph_helac_xsec_deltay_tot.Clone("graph_helac_xsec_deltay_tot")
    graph_helac_xsec_deltay_tot_cent.SetLineWidth(4)
    graph_helac_xsec_deltay_tot_cent.SetMarkerSize(0)
    for ibin in range(graph_helac_xsec_deltay_tot_cent.GetN()):
        graph_helac_xsec_deltay_tot_cent.SetPointError(ibin, graph_helac_xsec_deltay_tot_cent.GetErrorXlow(ibin),
                                                       graph_helac_xsec_deltay_tot_cent.GetErrorXhigh(ibin), 0., 0.)
    graph_helac_xsec_deltay_tot_cent.Draw("pz")
    graph_xsec_deltay_syst.Draw("2")
    graph_xsec_deltay_stat.Draw("pz")
    leg_pythia.SetY1NDC(0.38)
    leg_pythia.SetY2NDC(0.53)
    leg_pythia.SetX1NDC(0.60)
    leg_pythia.SetX2NDC(0.90)
    # leg_pythia.SetNColumns(1)
    leg_helac.SetY1NDC(0.56)
    leg_helac.SetY2NDC(0.74)
    leg_helac.SetX1NDC(0.60)
    leg_helac.SetX2NDC(0.90)
    leg_pythia.Draw()
    leg_helac.Draw()
    canv_xsec_deltay.SaveAs(
        f"dzero_jpsi_deltay_diff_xsection_pp13dot6TeV_HELACONIAsigmaeff{sigma_eff_helac:.0f}.pdf")

    canv_xsec_deltay_helaconly = ROOT.TCanvas("canv_xsec_deltay_helaconly", "", 1200, 1200)
    frame = canv_xsec_deltay_helaconly.DrawFrame(-2., 0., 5., 0.5, ";#Delta#it{y};d#sigma/d#Delta#it{y} (#mub)")
    frame.GetYaxis().SetDecimals()
    lat_label.DrawLatex(0.18, 0.9, "ALICE Preliminary")
    lat.DrawLatex(0.18, 0.86, "pp collisions, #sqrt{s} = 13.6 TeV, #it{L}_{int} = 39.7 pb^{#minus1}")
    lat.DrawLatex(0.18, 0.81, "Prompt D^{0}, |#it{y}| < 0.6, #it{p}_{T} > 0.5 GeV/#it{c}")
    lat.DrawLatex(0.18, 0.77, "Inclusive J/#Psi, 2.5 < #it{y} < 4.0, #it{p}_{T} > 0")
    graph_helac_xsec_deltay_sps.Draw("2")
    graph_helac_xsec_deltay_sps_cent.Draw("pz")
    graph_helac_xsec_deltay_dps.Draw("pz")
    graph_helac_xsec_deltay_tot.Draw("2")
    graph_helac_xsec_deltay_tot_cent.Draw("pz")
    graph_xsec_deltay_syst.Draw("2")
    graph_xsec_deltay_stat.Draw("pz")
    leg_helac.Draw()
    canv_xsec_deltay_helaconly.SaveAs(
        f"dzero_jpsi_deltay_diff_xsection_pp13dot6TeV_HELACONIAonly_HELACONIAsigmaeff{sigma_eff_helac:.0f}.pdf")

    canv_ratios = ROOT.TCanvas("canv_ratios", "", 1200, 600)
    canv_ratios_padassocprod = ROOT.TPad("canv_ratios_padassocprod", "canv_ratios_padassocprod", 0.64, 0., 1., 1.)
    canv_ratios_padassocprod.SetLeftMargin(0.)
    canv_ratios_padassocprod.Draw()
    canv_ratios_padassocprod.cd()
    frame_ratio_assocprod = ROOT.TH2D("frame_ratio_assocprod", ";#it{R}_{DJ/#Psi} (mb);", 1000, 1.2, 0.75e2, 3, 0.5, 3.5)
    frame_ratio_assocprod.GetYaxis().SetBinLabel(1, "")
    frame_ratio_assocprod.GetXaxis().SetMoreLogLabels()
    frame_ratio_assocprod.GetXaxis().SetLabelSize(0.058)
    frame_ratio_assocprod.GetXaxis().SetTitleSize(0.065)
    frame_ratio_assocprod.GetXaxis().SetLabelOffset(-0.01)
    frame_ratio_assocprod.GetXaxis().SetTitleOffset(0.75)
    frame_ratio_assocprod.GetXaxis().SetTickLength(0.018)
    frame_ratio_assocprod.GetYaxis().SetTickLength(0.05)
    frame_ratio_assocprod.GetXaxis().CenterTitle()
    frame_ratio_assocprod.DrawCopy()
    graph_pythia_ratio_tot.Draw("pz")
    graph_pythia_ratio_sps.Draw("pz")
    graph_pythia_ratio_dps.Draw("2")
    graph_ratio_syst.Draw("2")
    graph_ratio_stat.Draw("pz")
    leg_pythia.SetY1NDC(0.7)
    leg_pythia.SetY2NDC(0.92)
    leg_pythia.SetX1NDC(0.1)
    leg_pythia.SetX2NDC(0.82)
    leg_pythia.SetNColumns(2)
    leg_pythia.SetTextSize(0.06)
    leg_pythia.Draw()
    canv_ratios.cd()
    canv_ratios_padsinglehad = ROOT.TPad("canv_ratios_padsinglehad", "canv_ratios_padsinglehad", 0., 0., 0.64, 1.)
    canv_ratios_padsinglehad.SetLeftMargin(0.45)
    canv_ratios_padsinglehad.SetRightMargin(0.)
    canv_ratios_padsinglehad.Draw()
    canv_ratios_padsinglehad.cd()
    frame_ratio_singlehad = ROOT.TH2D("frame_ratio_singlehad", ";#sigma(J/#Psi)/#sigma(D^{0});",
                                      1000, 0.375e-3, 4.75e-2, 3, 0.5, 3.5)
    frame_ratio_singlehad.GetYaxis().SetBinLabel(
        2, "#splitline{ALICE, #sqrt{#it{s}} = 13.6 TeV}{#splitline{Prompt D^{0} |#it{y}|<0.6, #it{p}_{T}>0.5 GeV/#it{c}}{Inclusive J/#Psi 2.5<#it{y}<4, #it{p}_{T}>0}}") #
    frame_ratio_singlehad.GetYaxis().SetBinLabel(
        1, "#splitline{LHCb, #sqrt{#it{s}} = 7 TeV}{#splitline{Prompt D^{0} 2<#it{y}<4, 3<#it{p}_{T}<12 GeV/#it{c}}{Prompt J/#Psi 2<#it{y}<4, #it{p}_{T}<12 GeV/#it{c}}}") #
    frame_ratio_singlehad.GetYaxis().SetLabelSize(0.055)
    frame_ratio_singlehad.GetXaxis().SetLabelSize(0.042)
    frame_ratio_singlehad.GetXaxis().SetTitleSize(0.046)
    frame_ratio_singlehad.GetXaxis().SetTitleOffset(0.95)
    frame_ratio_singlehad.GetXaxis().SetMoreLogLabels()
    frame_ratio_singlehad.GetXaxis().CenterTitle()
    frame_ratio_singlehad.GetXaxis().SetNdivisions(508)
    frame_ratio_singlehad.DrawCopy()
    graph_pythia_ratio_singlehad.Draw("2")
    graph_ratio_singlehad_syst.Draw("2")
    graph_ratio_singlehad_stat.Draw("pz")
    lat_label.SetTextSize(0.05)
    lat.SetTextSize(0.045)
    lat_label.DrawLatex(0.48, 0.88, "ALICE Preliminary")
    lat.DrawLatex(0.48, 0.83, "pp collisions")
    lat.DrawLatex(0.48, 0.78, "#it{L}_{int} = 39.7 pb^{#minus1}")
    canv_ratios.SaveAs(f"dzero_jpsi_ratios_pp13dot6TeV_HELACONIAsigmaeff{sigma_eff_helac:.0f}.pdf")

    # evaluate chi2 with respect to HELAC-Onia
    chi2_exp, chi2_theor = 0., 0.
    ndf = graph_xsec_deltay_stat.GetN()
    for ibin in range(ndf):
        delta_y_cent = graph_xsec_deltay_stat.GetPointX(ibin)
        xsec_exp = graph_xsec_deltay_stat.GetPointY(ibin)
        xsec_stat = graph_xsec_deltay_stat.GetErrorYlow(ibin)
        xsec_syst = graph_xsec_deltay_syst.GetErrorYlow(ibin)
        for ibin_helac in range(graph_helac_xsec_deltay_tot.GetN()):
            if np.isclose(delta_y_cent, graph_helac_xsec_deltay_tot.GetPointX(ibin_helac)):
                xsec_helac_onia = graph_helac_xsec_deltay_tot.GetPointY(ibin_helac)
                if xsec_exp - xsec_helac_onia < 0:
                    xsec_helac_onia_unc = graph_helac_xsec_deltay_tot.GetErrorYlow(ibin_helac)
                else:
                    xsec_helac_onia_unc = graph_helac_xsec_deltay_tot.GetErrorYhigh(ibin_helac)
                break
        chi2_exp += (xsec_exp - xsec_helac_onia)**2 / (xsec_stat**2 + xsec_syst**2)
        chi2_theor += (xsec_exp - xsec_helac_onia)**2 / xsec_helac_onia_unc**2

    return chi2_exp / ndf, chi2_theor / ndf


def do_chi2_scan():
    """
    Method to extract sigma_eff via comparison with HELAC-Onia
    """

    sigma_eff_values = [10. + ivalue for ivalue in range(111)]

    graph_chi2_exp = ROOT.TGraph(1)
    graph_chi2_exp.SetName("graph_chi2_exp")
    graph_chi2_theory = ROOT.TGraph(1)
    graph_chi2_theory.SetName("graph_chi2_theory")
    set_style(graph_chi2_exp, ROOT.kAzure+4, ROOT.kFullCircle, 0.5, 0, 2, 1.)
    set_style(graph_chi2_theory, ROOT.kRed+1, ROOT.kFullSquare, 0.5, 0, 2, 1.)

    line_at_one = ROOT.TLine(10., 1., 120, 1.)
    line_at_one.SetLineWidth(1)
    line_at_one.SetLineColor(ROOT.kGray+1)
    line_at_one.SetLineStyle(9)

    for isigma_eff, sigma_eff in enumerate(sigma_eff_values):
        chi2_exp, chi2_theor = plot(sigma_eff)
        print(f"sigma_eff = {sigma_eff:.0f} mb, chi2(exp) = {chi2_exp:0.4f}, chi2(theory) = {chi2_theor:0.4f}")
        graph_chi2_exp.SetPoint(isigma_eff, sigma_eff, chi2_exp)
        graph_chi2_theory.SetPoint(isigma_eff, sigma_eff, chi2_theor)

    canv_chi2 = ROOT.TCanvas("canv_chi2", "", 500, 500)
    leg = ROOT.TLegend(0.45, 0.75, 0.8, 0.85)
    leg.SetTextSize(0.04)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetMargin(0.1)
    leg.AddEntry(graph_chi2_exp, "Experimental unc.", "lp")
    leg.AddEntry(graph_chi2_theory, "HELAC-Onia unc. (SPS)", "lp")
    canv_chi2.SetGridy()
    canv_chi2.SetGridx()
    frame = canv_chi2.DrawFrame(10., 0., 120., 5., ";#sigma_{eff} (mb);#chi^{2}")
    frame.GetYaxis().SetDecimals()
    line_at_one.Draw("same")
    graph_chi2_exp.Draw("lpz")
    graph_chi2_theory.Draw("lpz")
    leg.Draw()
    canv_chi2.SaveAs("chi2_sigmaeff_helac.pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("-s", "--sigma_eff_helac", type=float,
                        default=15., required=False,
                        help="sigma_eff parameter in mb in HELAC predictions")
    parser.add_argument("-e", "--extraction_sigma_eff", action="store_true",
                        default=False, required=False,
                        help="enable")
    args = parser.parse_args()

    if args.extraction_sigma_eff:
        do_chi2_scan()
    else:
        plot(args.sigma_eff_helac)
