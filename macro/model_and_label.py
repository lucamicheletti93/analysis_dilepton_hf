#this script include functions to recostruct the model and set the labels for the fitting part
import sys
import argparse
import yaml
import ROOT
import numpy as np
import os


def getLabel(val):
    """
    function to get a label for saving results without overwriting
    this is for just one value
    """
    label = ""
    if(val<0):
        label +="Neg"
    label += str(int(ROOT.TMath.Abs(val)))
    label += "p"
    label += str(int((10*ROOT.TMath.Abs(val))%10))
    return label

def getGlobalLabel(config):
    return config["inputs"]["inputLabel"]+("_weightedFit" if config['fit']['weighted'] else "")+"_dRap_min"+getLabel(config['fit']['min_dRap'])+"_max"+getLabel(config['fit']['max_dRap'])+"rapJpsi_min"+getLabel(config['fit']['min_jpsi_rap'])+"_max"+getLabel(config['fit']['max_jpsi_rap'])+"_ptJpsi_min"+getLabel(config['fit']['min_jpsi_pt'])+"_max"+getLabel(config['fit']['max_jpsi_pt'])+"_rapD0_min"+getLabel(config['fit']['min_d0_rap'])+"_max"+getLabel(config['fit']['max_d0_rap'])+"_ptD0_min"+getLabel(config['fit']['min_d0_pt'])+"_max"+getLabel(config['fit']['max_d0_pt'])
    
    
def constructWorkspace(config, includeJpsi, includeD0, nEvent, variations):
    '''
    This function reconstructs the model for all fits to avoid many copy-pastes in different functions
    '''
    if config["fit"]["JpsiChannel"] == "Jpsi2ee":
        titleSuffix = "e^{+}e^{-}"
    elif config["fit"]["JpsiChannel"] == "Jpsi2mumu":
        titleSuffix = "#mu^{+}#mu^{-}"
    else:
        print("Error: JpsiChannel not defined in the configuration file.")
        sys.exit(1)

    wsLabel = "ws_"
    if includeJpsi: wsLabel += "Jpsi"
    if includeD0:  wsLabel += "D0"
    wsLabel += getGlobalLabel(config)

    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
    ws = ROOT.RooWorkspace(wsLabel)

    minFitRangeD0 = config["fit"]["min_fit_range_d0"]
    maxFitRangeD0 = config["fit"]["max_fit_range_d0"]
    minFitRangeJpsi = config["fit"]["min_fit_range_jpsi"]
    maxFitRangeJpsi = config["fit"]["max_fit_range_jpsi"]

    minFitRangeD0 = minFitRangeD0 + config["fit"]["var_min_fit_range_d0"][variations[0]]
    maxFitRangeD0 = maxFitRangeD0 + config["fit"]["var_max_fit_range_d0"][variations[0]]
    minFitRangeJpsi = minFitRangeJpsi + config["fit"]["var_min_fit_range_jpsi"][variations[1]]
    maxFitRangeJpsi = maxFitRangeJpsi + config["fit"]["var_max_fit_range_jpsi"][variations[1]]

    print(f'[INFO] Fit Range D0: {minFitRangeD0} - {maxFitRangeD0}')
    print(f'[INFO] Fit Range Jpsi: {minFitRangeJpsi} - {maxFitRangeJpsi}')

    mD0   = ROOT.RooRealVar("fMassDmes", "#it{m}_{#piK} (GeV/#it{c}^{2})", minFitRangeD0, maxFitRangeD0); getattr(ws, "import")(mD0)
    mJpsi = ROOT.RooRealVar("fMass", f"#it{{m}}_{{{titleSuffix}}} (GeV/#it{{c}}^{{2}})", minFitRangeJpsi, maxFitRangeJpsi); getattr(ws, "import")(mJpsi)
    ptD0   = ROOT.RooRealVar("fPtDmes", "#it{p}_{T,#piK} (GeV/#it{c}^{2})", config["fit"]["min_d0_pt"], config["fit"]["max_d0_pt"]); getattr(ws, "import")(ptD0)
    ptJpsi = ROOT.RooRealVar("fPtJpsi", f"#it{{p}}_{{T,{titleSuffix}}} (GeV/#it{{c}}^{{2}})", config["fit"]["min_jpsi_pt"], config["fit"]["max_jpsi_pt"]); getattr(ws, "import")(ptJpsi)
    dRap = ROOT.RooRealVar("fDeltaY", f"y_{{{titleSuffix}}} - y_{{#piK}}", config["fit"]["min_dRap"], config["fit"]["max_dRap"]); getattr(ws, "import")(dRap)
    rapJpsi = ROOT.RooRealVar("fRapJpsi", f"y_{{{titleSuffix}}}", config["fit"]["min_jpsi_rap"], config["fit"]["max_jpsi_rap"]); getattr(ws, "import")(rapJpsi)
    if config["fit"]["weighted"]:
        weight = ROOT.RooRealVar("weight", "weight", 0, 1e10); getattr(ws, "import")(weight)
            
    # Pdfs
    meanJpsi  = ROOT.RooRealVar(config["fit"]["cb_par_jpsi_name"][0], config["fit"]["cb_par_jpsi_name"][0], config["fit"]["cb_par_jpsi_val"][0], config["fit"]["cb_par_jpsi_lw_lim"][0], config["fit"]["cb_par_jpsi_up_lim"][0]); meanJpsi.setConstant(config["fit"]["cb_par_jpsi_is_const"][0])
    sigmaJpsi = ROOT.RooRealVar(config["fit"]["cb_par_jpsi_name"][1], config["fit"]["cb_par_jpsi_name"][1], config["fit"]["cb_par_jpsi_val"][1], config["fit"]["cb_par_jpsi_lw_lim"][1], config["fit"]["cb_par_jpsi_up_lim"][1]); sigmaJpsi.setConstant(config["fit"]["cb_par_jpsi_is_const"][1])
    alphaJpsi = ROOT.RooRealVar(config["fit"]["cb_par_jpsi_name"][2], config["fit"]["cb_par_jpsi_name"][2], config["fit"]["cb_par_jpsi_val"][2], config["fit"]["cb_par_jpsi_lw_lim"][2], config["fit"]["cb_par_jpsi_up_lim"][2]); alphaJpsi.setConstant(config["fit"]["cb_par_jpsi_is_const"][2])
    nJpsi     = ROOT.RooRealVar(config["fit"]["cb_par_jpsi_name"][3], config["fit"]["cb_par_jpsi_name"][3], config["fit"]["cb_par_jpsi_val"][3], config["fit"]["cb_par_jpsi_lw_lim"][3], config["fit"]["cb_par_jpsi_up_lim"][3]); nJpsi.setConstant(config["fit"]["cb_par_jpsi_is_const"][3])
    cbPdfJpsi = ROOT.RooCBShape("cbPdfJpsi", "Crystal Ball J/psi", mJpsi, meanJpsi, sigmaJpsi, alphaJpsi, nJpsi)

    ## adding the psi2S
    if config["fit"]["add_psi2s"]:
        meanPsi2S  = ROOT.RooRealVar(config["fit"]["cb_par_psi2s_name"][0], config["fit"]["cb_par_psi2s_name"][0], config["fit"]["cb_par_psi2s_val"][0], config["fit"]["cb_par_psi2s_lw_lim"][0], config["fit"]["cb_par_psi2s_up_lim"][0]); meanPsi2S.setConstant(config["fit"]["cb_par_psi2s_is_const"][0])
        sigmaPsi2S = ROOT.RooRealVar(config["fit"]["cb_par_psi2s_name"][1], config["fit"]["cb_par_psi2s_name"][1], config["fit"]["cb_par_psi2s_val"][1], config["fit"]["cb_par_psi2s_lw_lim"][1], config["fit"]["cb_par_psi2s_up_lim"][1]); sigmaPsi2S.setConstant(config["fit"]["cb_par_psi2s_is_const"][1])
        alphaPsi2S = ROOT.RooRealVar(config["fit"]["cb_par_psi2s_name"][2], config["fit"]["cb_par_psi2s_name"][2], config["fit"]["cb_par_psi2s_val"][2], config["fit"]["cb_par_psi2s_lw_lim"][2], config["fit"]["cb_par_psi2s_up_lim"][2]); alphaPsi2S.setConstant(config["fit"]["cb_par_psi2s_is_const"][2])
        nPsi2S     = ROOT.RooRealVar(config["fit"]["cb_par_psi2s_name"][3], config["fit"]["cb_par_psi2s_name"][3], config["fit"]["cb_par_psi2s_val"][3], config["fit"]["cb_par_psi2s_lw_lim"][3], config["fit"]["cb_par_psi2s_up_lim"][3]); nPsi2S.setConstant(config["fit"]["cb_par_psi2s_is_const"][3])
        cbPdfPsi2S = ROOT.RooCBShape("cbPdfPsi2S", "Crystal Ball J/psi", mJpsi, meanPsi2S, sigmaJpsi, alphaPsi2S, nPsi2S)

    ##same bkg for both
    chebyParsJpsi = [ROOT.RooRealVar(config["fit"]["cheby_par_jpsi_name"][i], config["fit"]["cheby_par_jpsi_name"][i], config["fit"]["cheby_par_jpsi_val"][i], config["fit"]["cheby_par_jpsi_lw_lim"][i], config["fit"]["cheby_par_jpsi_up_lim"][i]) for i in range(3)]
    for i in range(0, 3): chebyParsJpsi[i].setConstant(config["fit"]["cheby_par_jpsi_is_const"][i])
    chebyPdfJpsi  = ROOT.RooChebychev("chebyPdfJpsi", "Cheby for Bkg1", mJpsi, ROOT.RooArgList(*chebyParsJpsi))

    meanD0  = ROOT.RooRealVar(config["fit"]["cb_par_d0_name"][0], config["fit"]["cb_par_d0_name"][0], config["fit"]["cb_par_d0_val"][0], config["fit"]["cb_par_d0_lw_lim"][0], config["fit"]["cb_par_d0_up_lim"][0]); meanD0.setConstant(config["fit"]["cb_par_d0_is_const"][0])
    sigmaD0 = ROOT.RooRealVar(config["fit"]["cb_par_d0_name"][1], config["fit"]["cb_par_d0_name"][1], config["fit"]["cb_par_d0_val"][1], config["fit"]["cb_par_d0_lw_lim"][1], config["fit"]["cb_par_d0_up_lim"][1]); sigmaD0.setConstant(config["fit"]["cb_par_d0_is_const"][1])
    alphaD0 = ROOT.RooRealVar(config["fit"]["cb_par_d0_name"][2], config["fit"]["cb_par_d0_name"][2], config["fit"]["cb_par_d0_val"][2], config["fit"]["cb_par_d0_lw_lim"][2], config["fit"]["cb_par_d0_up_lim"][1]); alphaD0.setConstant(config["fit"]["cb_par_d0_is_const"][2])
    nD0     = ROOT.RooRealVar(config["fit"]["cb_par_d0_name"][3], config["fit"]["cb_par_d0_name"][3], config["fit"]["cb_par_d0_val"][3], config["fit"]["cb_par_d0_lw_lim"][3], config["fit"]["cb_par_d0_up_lim"][2]); nD0.setConstant(config["fit"]["cb_par_d0_is_const"][3])
    cbPdfD0 = ROOT.RooCBShape("cbPdfD0", "Crystal Ball D0", mD0, meanD0, sigmaD0, alphaD0, nD0)

    chebyParsD0 = [ROOT.RooRealVar(config["fit"]["cheby_par_d0_name"][i], config["fit"]["cheby_par_d0_name"][i], config["fit"]["cheby_par_d0_val"][i], config["fit"]["cheby_par_d0_lw_lim"][i], config["fit"]["cheby_par_d0_up_lim"][i]) for i in range(3)]
    for i in range(0, 3): chebyParsD0[i].setConstant(config["fit"]["cheby_par_d0_is_const"][i])
    chebyPdfD0  = ROOT.RooChebychev("chebyPdfD0", "Cheby for Bkg1", mD0, ROOT.RooArgList(*chebyParsD0))

    # reflection
    fRefl = ROOT.TFile(config["reflections"]["data"], "READ")
    hRefl = fRefl.Get(config["reflections"]["refl"])
    hSig = fRefl.Get(config["reflections"]["signal"])
    bin_min_refl = hRefl.FindBin(minFitRangeD0)
    bin_max_refl = hRefl.FindBin(maxFitRangeD0)
    bin_min_sig = hSig.FindBin(minFitRangeD0)
    bin_max_sig = hSig.FindBin(maxFitRangeD0)

    integral_hRefl = hRefl.Integral(bin_min_refl, bin_max_refl)
    integral_hSig = hSig.Integral(bin_min_sig, bin_max_sig)

    if integral_hSig != 0:
        reflRatio = integral_hRefl / integral_hSig
    else:
        print(f"Error: The integral of hSig in the mass range is zero.")
        sys.exit(1)  # Exit the script with a non-zero status to indicate an error

    reflFracD0 = ROOT.RooRealVar("reflFracD0", "reflFracD0", reflRatio); reflFracD0.setConstant(True)
    
    refl_template_D0 = ROOT.RooDataHist("refl_template_D0", "refl_template_D0", ROOT.RooArgList(mD0), hRefl)
    templatePdfReflD0 = ROOT.RooHistPdf("templatePdfReflD0", "templatePdfReflD0", ROOT.RooArgList(mD0), refl_template_D0)

    
#### now make the products of the pdfs depending on the case
    if (includeJpsi and not includeD0):
        nSigJpsi  = ROOT.RooRealVar("nSigJpsi", "Jpsi signal", 0.03*nEvent,  0.0001*nEvent, nEvent) ## to be changed to not be hardcoded
        nSigPsi2S  = ROOT.RooRealVar("nSigPsi2S", "psi(2S) signal", 0.00003*nEvent, 0.000001*nEvent, 0.0003*nEvent) ## to be changed to not be hardcoded
        nBkgJpsi  = ROOT.RooRealVar("nBkgJpsi", "Jpsi background", 0.07*nEvent, 0.0001*nEvent, nEvent) ## to be changed to not be hardcoded
        if config["fit"]["add_psi2s"]:
            model = ROOT.RooAddPdf("model", "sigJpsi + nSigPsi2S + bkgJpsi", ROOT.RooArgList(cbPdfJpsi, cbPdfPsi2S, chebyPdfJpsi), ROOT.RooArgList(nSigJpsi, nSigPsi2S, nBkgJpsi))
        else:
            model = ROOT.RooAddPdf("model", "sigJpsi + bkgJpsi", ROOT.RooArgList(cbPdfJpsi, chebyPdfJpsi), ROOT.RooArgList(nSigJpsi, nBkgJpsi))
    elif (includeD0 and not includeJpsi):
        nSigD0  = ROOT.RooRealVar("nSigD0", "D0 signal", 0.003*nEvent, 0.0003*nEvent, 0.5*nEvent) ## to be changed to not be hardcoded
        nBkgD0  = ROOT.RooRealVar("nBkgD0", "D0 background", 0.05*nEvent, 0.003*nEvent, 0.5*nEvent) ## to be changed to not be hardcoded
        nReflD0  = ROOT.RooFormulaVar("nReflD0", "reflFracD0 * nSigD0", ROOT.RooArgList(reflFracD0,nSigD0))
        model = ROOT.RooAddPdf("model", "sigD0 + bkgD0 + reflD0", ROOT.RooArgList(cbPdfD0, chebyPdfD0, templatePdfReflD0), ROOT.RooArgList(nSigD0, nBkgD0, nReflD0))
    elif(includeD0 and includeJpsi):
        nJPsiD0  = ROOT.RooRealVar("nJPsiD0", "number of JPsi-D0", config["fit"]["norm_par_sig_val"][0], config["fit"]["norm_par_sig_lw_lim"][0], config["fit"]["norm_par_sig_up_lim"][0])
        nBkgJPsi = ROOT.RooRealVar("nBkgJPsi", "number of JPsi-Bkg", config["fit"]["norm_par_sig_val"][1], config["fit"]["norm_par_sig_lw_lim"][1], config["fit"]["norm_par_sig_up_lim"][1])
        nBkgD0   = ROOT.RooRealVar("nBkgD0", "number of D0-Bkg", config["fit"]["norm_par_sig_val"][2], config["fit"]["norm_par_sig_lw_lim"][2], config["fit"]["norm_par_sig_up_lim"][2])
        nBkgBkg  = ROOT.RooRealVar("nBkgBkg", "number of Bkg-Bkg", config["fit"]["norm_par_sig_val"][3], config["fit"]["norm_par_sig_lw_lim"][3], config["fit"]["norm_par_sig_up_lim"][3])
        if config["fit"]["add_psi2s"]:
            rPsi2SJPsi  = ROOT.RooRealVar("rPsi2SJPsi", "Psi2S/JPsi ratio", config["fit"]["norm_par_sig_val"][4], config["fit"]["norm_par_sig_lw_lim"][4], config["fit"]["norm_par_sig_up_lim"][4])
            rPsi2SJPsi.setConstant(config["fit"]["norm_par_sig_is_const"][4])
            nPsi2SD0 = ROOT.RooFormulaVar("nPsi2SD0", "rPsi2SJPsi * nJPsiD0", ROOT.RooArgList(rPsi2SJPsi,nJPsiD0))
            nBkgPsi2S = ROOT.RooFormulaVar("nBkgPsi2S", "rPsi2SJPsi * nBkgJPsi", ROOT.RooArgList(rPsi2SJPsi,nBkgJPsi))

        nReflJPsi = ROOT.RooFormulaVar("nReflJPsi", "reflFracD0 * nJPsiD0", ROOT.RooArgList(reflFracD0,nJPsiD0))
        if config["fit"]["add_psi2s"]:
            nReflPsi2S = ROOT.RooFormulaVar("nReflPsi2S", "reflFracD0 * nPsi2SD0", ROOT.RooArgList(reflFracD0,nPsi2SD0))
        nReflBkg  = ROOT.RooFormulaVar("nReflBkg", "reflFracD0 * nBkgD0", ROOT.RooArgList(reflFracD0,nBkgD0))
    
        # Product of Pdfs
        sigD0_sigJpsiPdf = ROOT.RooProdPdf("sigD0_sigJpsiPdf", "Jpsi Pdf * D0 Pdf", ROOT.RooArgList(cbPdfD0, cbPdfJpsi))
        bkgD0_sigJpsiPdf = ROOT.RooProdPdf("bkgD0_sigJpsiPdf", "Jpsi Pdf * Bkg2 Pdf", ROOT.RooArgList(chebyPdfD0, cbPdfJpsi))
        reflD0_sigJpsiPdf = ROOT.RooProdPdf("reflD0_sigJpsiPdf", "Jpsi Pdf * refl Pdf", ROOT.RooArgList(templatePdfReflD0, cbPdfJpsi))
        if config["fit"]["add_psi2s"]:
            sigD0_sigPsi2SPdf = ROOT.RooProdPdf("sigD0_sigPsi2SPdf", "Psi2S Pdf * D0 Pdf", ROOT.RooArgList(cbPdfD0, cbPdfPsi2S))
            bkgD0_sigPsi2SPdf = ROOT.RooProdPdf("bkgD0_sigPsi2SPdf", "Psi2S Pdf * Bkg2 Pdf", ROOT.RooArgList(chebyPdfD0, cbPdfPsi2S))
            reflD0_sigPsi2SPdf = ROOT.RooProdPdf("reflD0_sigPsi2SPdf", "Psi2S Pdf * refl Pdf", ROOT.RooArgList(templatePdfReflD0, cbPdfPsi2S))
    
        bkgJpsi_sigD0Pdf = ROOT.RooProdPdf("bkgJpsi_sigD0Pdf", "D0 Pdf * Bkg1 Pdf", ROOT.RooArgList(chebyPdfJpsi, cbPdfD0))
        bkgD0_bkgJpsiPdf = ROOT.RooProdPdf("bkgD0_bkgJpsiPdf", "Bkg1 Pdf * Bkg2 Pdf", ROOT.RooArgList(chebyPdfD0, chebyPdfJpsi))
        reflD0_bkgJpsiPdf = ROOT.RooProdPdf("reflD0_bkgJpsiPdf", "Bkg1 Pdf * refl Pdf", ROOT.RooArgList(templatePdfReflD0, chebyPdfJpsi))
        
        # Sum all Pdfs
        if config["fit"]["add_psi2s"]:
            model = ROOT.RooAddPdf("model", "sigD0_sigJpsiPdf + bkgD0_sigJpsiPdf + reflD0_sigJpsiPdf + sigD0_sigPsi2SPdf + bkgD0_sigPsi2SPdf + reflD0_sigPsi2SPdf + bkgJpsi_sigD0Pdf + bkgD0_bkgJpsiPdf + reflD0_bkgJpsiPdf", ROOT.RooArgList(sigD0_sigJpsiPdf, bkgD0_sigJpsiPdf, reflD0_sigJpsiPdf, sigD0_sigPsi2SPdf, bkgD0_sigPsi2SPdf, reflD0_sigPsi2SPdf,bkgJpsi_sigD0Pdf, bkgD0_bkgJpsiPdf, reflD0_bkgJpsiPdf), ROOT.RooArgList(nJPsiD0, nBkgJPsi, nReflJPsi, nPsi2SD0, nBkgPsi2S, nReflPsi2S,nBkgD0, nBkgBkg, nReflBkg))
        else:
            model = ROOT.RooAddPdf("model", "sigD0_sigJpsiPdf + bkgD0_sigJpsiPdf + reflD0_sigJpsiPdf  + bkgJpsi_sigD0Pdf + bkgD0_bkgJpsiPdf + reflD0_bkgJpsiPdf", ROOT.RooArgList(sigD0_sigJpsiPdf, bkgD0_sigJpsiPdf, reflD0_sigJpsiPdf,bkgJpsi_sigD0Pdf, bkgD0_bkgJpsiPdf, reflD0_bkgJpsiPdf), ROOT.RooArgList(nJPsiD0, nBkgJPsi, nReflJPsi,nBkgD0, nBkgBkg, nReflBkg))

    else:
        print("Error: you need to include at least one of the two dimensions")
        return 0

    getattr(ws, "import")(model)
    return ws

