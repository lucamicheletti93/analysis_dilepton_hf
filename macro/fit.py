import sys
import argparse
import yaml
import ROOT
sys.path.append('../utils')
from plot_library import LoadStyle, SetGraStat, SetGraSyst, SetLegend

def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='config.yml', help='config file name')
    parser.add_argument("--do_prefit", help="Do the pre-fit to data", action="store_true")
    parser.add_argument("--do_fit", help="Do the fit to data / toy MC", action="store_true")
    parser.add_argument("--do_upper_limit", help="Do the calculation of the upper limit", action="store_true")
    parser.add_argument("--do_prefilter", help="Apply selections", action="store_true")
    parser.add_argument("--do_plot_results", help="Plot results", action="store_true")
    args = parser.parse_args()
    

    with open(args.cfgFileName, 'r') as yml_cfg:
        inputCfg = yaml.load(yml_cfg, yaml.FullLoader)

    if args.do_prefit:
        prefit()

    if args.do_fit:
        fit(inputCfg)

    if args.do_upper_limit:
        upper_limit()

    if args.do_prefilter:
        prefilter(inputCfg)

    if args.do_plot_results:
        plot_results()

def prefilter(config):
    """
    function to apply pre-filters on D0 tree according to BDT output (see config_fit.yml)
    """
    fIn = ROOT.TFile(config["prefilter"]["data"], "READ")
    tree = fIn.Get(config["prefilter"]["tree"])
    rDataFrame = ROOT.RDataFrame(tree).Filter(config["prefilter"]["cuts"])
    fOutName = config["prefilter"]["data"].replace(".root", config["prefilter"]["suffix"] + ".root")
    rDataFrame.Snapshot(config["prefilter"]["tree"], fOutName)
    fIn.Close()

def prefit():
    """
    function to fit the integrated spectrum for D0 and J/psi to extract shape parameters
    """
    LoadStyle()
    ROOT.gStyle.SetPalette(ROOT.kBird)

    letexTitle = ROOT.TLatex()
    letexTitle.SetTextSize(0.040)
    letexTitle.SetNDC()
    letexTitle.SetTextFont(42)

    # Variables
    mD0   = ROOT.RooRealVar("D-meson mass", "#it{m}_{#pi#it{K}} (GeV/#it{c}^{2})", 1.76, 2.00)
    mJpsi = ROOT.RooRealVar("Dimuon mass", "#it{m}_{#mu#mu} (GeV/#it{c}^{2})", 2.50, 4.00)

    # Pre-fit D0
    # Load data for pre-fit for D0 - taken from BDT selection applied to TH3D
    fInD0 = ROOT.TFile("/Users/lucamicheletti/cernbox/Jpsi_D0_analysis/train_137585/massD0.root", "READ")    
    hMassDmeson = fInD0.Get("hist_mass_d0_pt0_50")
    dataHistD0 = ROOT.RooDataHist("dataHistD0", "dataHistD0", [mD0], Import=hMassDmeson)
    
    meanD0  = ROOT.RooRealVar("meanD0", "meanD0", 1.850, 1.750, 1.950)
    sigmaD0 = ROOT.RooRealVar("sigmaD0", "sigmaD0", 0.050, 0.010, 0.080)
    alphaD0 = ROOT.RooRealVar("alphaD0", "alphaD0", 0.500, 0.000, 5.000)
    nD0     = ROOT.RooRealVar("nD0", "nD0", 1.000, 0.000, 10.000)
    cbPdfD0 = ROOT.RooCBShape("cbPdfD0", "Crystal Ball J/psi", mD0, meanD0, sigmaD0, alphaD0, nD0)

    chebyParsD0 = [ROOT.RooRealVar(f"cheb_coeff_{i}", f"Coeff_{i}", 0.05, -10, 10) for i in range(3)]
    chebyPdfD0 = ROOT.RooChebychev("chebyPdfD0", "Cheby for Bkg1", mD0, ROOT.RooArgList(*chebyParsD0))

    nSigD0  = ROOT.RooRealVar("nSigD0", "D0 signal", 1e6, 0., 1e10)
    nBkgD0  = ROOT.RooRealVar("nBkgD0", "D0 background", 1e6, 0., 1e10)
    modelD0 = ROOT.RooAddPdf("modelD0", "sigD0 + bkgD0", ROOT.RooArgList(cbPdfD0, chebyPdfD0), ROOT.RooArgList(nSigD0, nBkgD0))

    fitResultD0 = modelD0.fitTo(dataHistD0, ROOT.RooFit.PrintLevel(3), ROOT.RooFit.Optimize(1), ROOT.RooFit.Hesse(1), ROOT.RooFit.Minos(1), ROOT.RooFit.Strategy(2), ROOT.RooFit.Save(1))

    mD0frame = mD0.frame(Title=" ")
    dataHistD0.plotOn(mD0frame)
    modelD0.plotOn(mD0frame)
    modelD0.plotOn(mD0frame, Name={"Sig"}, Components={cbPdfD0}, LineStyle="--", LineColor=ROOT.kRed+1)
    modelD0.plotOn(mD0frame, Name={"Bkg"}, Components={chebyPdfD0}, LineStyle="--", LineColor=ROOT.kAzure+4)

    # Pre-fit J/psi
    # Load data for pre-fit for J/psi - Taken from AnalysisResults
    fInJpsi = ROOT.TFile("/Users/lucamicheletti/cernbox/Jpsi_D0_analysis/train_137585/AnalysisResults_LHC22_pass4_highIR.root", "READ")
    hlistInJpsi = fInJpsi.Get("table-maker-jpsi-hf/output")
    listInJpsi = hlistInJpsi.FindObject("JPsi")
    hMassJpsiTmp = listInJpsi.FindObject("hMassVsPtJPsi")
    hMassJPsi = hMassJpsiTmp.ProjectionY("hMassJPsi")

    dataHistJpsi = ROOT.RooDataHist("dataHistJpsi", "dataHistJpsi", [mJpsi], Import=hMassJPsi)

    meanJpsi  = ROOT.RooRealVar("meanJpsi", "meanJpsi", 3.096, 2.900, 3.250)
    sigmaJpsi = ROOT.RooRealVar("sigmaJpsi", "sigmaJpsi", 0.090, 0.060, 0.120)
    alphaJpsi = ROOT.RooRealVar("alphaJpsi", "alphaJpsi", 0.500, 0.000, 5.000)
    nJpsi     = ROOT.RooRealVar("nJpsi", "nJpsi", 1.000, 0.000, 10.000)
    cbPdfJpsi = ROOT.RooCBShape("cbPdfJpsi", "Crystal Ball J/psi", mJpsi, meanJpsi, sigmaJpsi, alphaJpsi, nJpsi)

    chebyParsJpsi = [ROOT.RooRealVar(f"cheb_coeff_{i}", f"Coeff_{i}", 0.00, -100, 100) for i in range(3)]
    chebyPdfJpsi = ROOT.RooChebychev("chebyPdfJpsi", "Cheby for Bkg1", mJpsi, ROOT.RooArgList(*chebyParsJpsi))

    nSigJpsi  = ROOT.RooRealVar("nSigJpsi", "Jpsi signal", 1e6, 0., 1e8)
    nBkgJpsi  = ROOT.RooRealVar("nBkgJpsi", "Jpsi background", 1e6, 0., 1e8)
    modelJpsi = ROOT.RooAddPdf("modelJpsi", "sigJpsi + bkgJpsi", ROOT.RooArgList(cbPdfJpsi, chebyPdfJpsi), ROOT.RooArgList(nSigJpsi, nBkgJpsi))

    fitResultJpsi = modelJpsi.fitTo(dataHistJpsi, ROOT.RooFit.PrintLevel(3), ROOT.RooFit.Optimize(1), ROOT.RooFit.Hesse(1), ROOT.RooFit.Minos(1), ROOT.RooFit.Strategy(2), ROOT.RooFit.Save(1))

    mJpsiframe = mJpsi.frame(Title=" ")
    dataHistJpsi.plotOn(mJpsiframe)
    modelJpsi.plotOn(mJpsiframe)
    modelJpsi.plotOn(mJpsiframe, Name={"Sig"}, Components={cbPdfJpsi}, LineStyle="--", LineColor=ROOT.kRed+1)
    modelJpsi.plotOn(mJpsiframe, Name={"Bkg"}, Components={chebyPdfJpsi}, LineStyle="--", LineColor=ROOT.kAzure+4)

    canvasFitD0 = ROOT.TCanvas("canvasFitD0", "canvasFitD0", 800, 600)
    canvasFitD0.SetTickx(1)
    canvasFitD0.SetTicky(1)
    mD0frame.GetYaxis().SetTitleOffset(1.4)
    mD0frame.Draw()

    letexTitle.DrawLatex(0.65, 0.74, "#it{N}_{D0} = %1.0f #pm %1.0f" % (nSigD0.getVal(), nSigD0.getError()))
    letexTitle.DrawLatex(0.65, 0.68, "#it{#mu}_{D0} = %4.3f #pm %4.3f" % (meanD0.getVal(), meanD0.getError()))
    letexTitle.DrawLatex(0.65, 0.62, "#it{#sigma}_{D0} = %4.3f #pm %4.3f" % (sigmaD0.getVal(), sigmaD0.getError()))
    
    canvasFitD0.Update()

    canvasFitJpsi = ROOT.TCanvas("canvasFitJpsi", "canvasFitJpsi", 800, 600)
    canvasFitJpsi.SetTickx(1)
    canvasFitJpsi.SetTicky(1)
    mJpsiframe.GetYaxis().SetTitleOffset(1.4)
    mJpsiframe.Draw()

    letexTitle.DrawLatex(0.65, 0.74, "#it{N}_{J/#psi} = %1.0f #pm %1.0f" % (nSigJpsi.getVal(), nSigJpsi.getError()))
    letexTitle.DrawLatex(0.65, 0.68, "#it{#mu}_{J/#psi} = %4.3f #pm %4.3f" % (meanJpsi.getVal(), meanJpsi.getError()))
    letexTitle.DrawLatex(0.65, 0.62, "#it{#sigma}_{J/#psi} = %4.3f #pm %4.3f" % (sigmaJpsi.getVal(), sigmaJpsi.getError()))
    
    canvasFitJpsi.Update()

    canvasFitD0.SaveAs("figures/prefitD0_new.pdf")
    canvasFitJpsi.SaveAs("figures/prefitJpsi_new.pdf")

    fitResultD0.Print()
    fitResultJpsi.Print()
    input()
    exit()

def fit(config):
    """
    function to fit the 2D distribution of D0 - J/psi masses
    - toy_mc:   True/False (enable the generation of toy MC sample to validate the result)
    - unbinned: True/False (enable the fit of unbinned datasets)
    """
    # Variables
    minFitRangeD0 = config["fit"]["min_fit_range_d0"]
    maxFitRangeD0 = config["fit"]["max_fit_range_d0"]
    minFitRangeJpsi = config["fit"]["min_fit_range_jpsi"]
    maxFitRangeJpsi = config["fit"]["max_fit_range_jpsi"]

    mD0   = ROOT.RooRealVar("fMassD0", "#it{m}_{#pi#it{K}} (GeV/#it{c}^{2})", minFitRangeD0, maxFitRangeD0)
    mJpsi = ROOT.RooRealVar("fMass", "#it{m}_{#mu#mu} (GeV/#it{c}^{2})", minFitRangeJpsi, maxFitRangeJpsi)

    # Yields parameters
    genJpsiD0  = config["fit"]["norm_par_sig_val"][0]
    genPsi2SD0  = config["fit"]["norm_par_sig_val"][1]
    genBkgJpsi = config["fit"]["norm_par_sig_val"][2]
    genBkgPsi2S = config["fit"]["norm_par_sig_val"][3]
    genBkgD0   = config["fit"]["norm_par_sig_val"][4]
    genBkgBkg  = config["fit"]["norm_par_sig_val"][5]
    

    nJPsiD0  = ROOT.RooRealVar("nJPsiD0", "number of JPsi-D0", genJpsiD0, config["fit"]["norm_par_sig_lw_lim"][0], config["fit"]["norm_par_sig_up_lim"][0])
    nPsi2SD0  = ROOT.RooRealVar("nPsi2SD0", "number of Psi2S-D0", genPsi2SD0, config["fit"]["norm_par_sig_lw_lim"][1], config["fit"]["norm_par_sig_up_lim"][1])
    nBkgJPsi = ROOT.RooRealVar("nBkgJPsi", "number of JPsi-Bkg", genBkgJpsi, config["fit"]["norm_par_sig_lw_lim"][2], config["fit"]["norm_par_sig_up_lim"][2])
    nBkgPsi2S = ROOT.RooRealVar("nBkgPsi2S", "number of Psi2S-Bkg", genBkgPsi2S, config["fit"]["norm_par_sig_lw_lim"][3], config["fit"]["norm_par_sig_up_lim"][3])
    nBkgD0   = ROOT.RooRealVar("nBkgD0", "number of D0-Bkg", genBkgD0, config["fit"]["norm_par_sig_lw_lim"][4], config["fit"]["norm_par_sig_up_lim"][4])
    nBkgBkg  = ROOT.RooRealVar("nBkgBkg", "number of Bkg-Bkg", genBkgBkg, config["fit"]["norm_par_sig_lw_lim"][5], config["fit"]["norm_par_sig_up_lim"][5])
    
    # Pdfs
    meanJpsi  = ROOT.RooRealVar(config["fit"]["cb_par_jpsi_name"][0], config["fit"]["cb_par_jpsi_name"][0], config["fit"]["cb_par_jpsi_val"][0], config["fit"]["cb_par_jpsi_lw_lim"][0], config["fit"]["cb_par_jpsi_up_lim"][0]); meanJpsi.setConstant(config["fit"]["cb_par_jpsi_is_const"][0])
    sigmaJpsi = ROOT.RooRealVar(config["fit"]["cb_par_jpsi_name"][1], config["fit"]["cb_par_jpsi_name"][1], config["fit"]["cb_par_jpsi_val"][1], config["fit"]["cb_par_jpsi_lw_lim"][1], config["fit"]["cb_par_jpsi_up_lim"][1]); sigmaJpsi.setConstant(config["fit"]["cb_par_jpsi_is_const"][1])
    alphaJpsi = ROOT.RooRealVar(config["fit"]["cb_par_jpsi_name"][2], config["fit"]["cb_par_jpsi_name"][2], config["fit"]["cb_par_jpsi_val"][2], config["fit"]["cb_par_jpsi_lw_lim"][2], config["fit"]["cb_par_jpsi_up_lim"][2]); alphaJpsi.setConstant(config["fit"]["cb_par_jpsi_is_const"][2])
    nJpsi     = ROOT.RooRealVar(config["fit"]["cb_par_jpsi_name"][3], config["fit"]["cb_par_jpsi_name"][3], config["fit"]["cb_par_jpsi_val"][3], config["fit"]["cb_par_jpsi_lw_lim"][3], config["fit"]["cb_par_jpsi_up_lim"][3]); nJpsi.setConstant(config["fit"]["cb_par_jpsi_is_const"][3])
    cbPdfJpsi = ROOT.RooCBShape("cbPdfJpsi", "Crystal Ball J/psi", mJpsi, meanJpsi, sigmaJpsi, alphaJpsi, nJpsi)

    ## adding the psi2S
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
    bin_min_refl = hRefl.FindBin(config["fit"]["min_fit_range_d0"])
    bin_max_refl = hRefl.FindBin(config["fit"]["max_fit_range_d0"])
    bin_min_sig = hSig.FindBin(config["fit"]["min_fit_range_d0"])
    bin_max_sig = hSig.FindBin(config["fit"]["max_fit_range_d0"])

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
    nReflJPsi = ROOT.RooFormulaVar("nReflJPsi", "reflFracD0 * nJPsiD0", ROOT.RooArgList(reflFracD0,nJPsiD0))
    nReflPsi2S = ROOT.RooFormulaVar("nReflPsi2S", "reflFracD0 * nPsi2SD0", ROOT.RooArgList(reflFracD0,nPsi2SD0))
    nReflBkg  = ROOT.RooFormulaVar("nReflBkg", "reflFracD0 * nBkgD0", ROOT.RooArgList(reflFracD0,nBkgD0))
    
    # Product of Pdfs
    sigD0_sigJpsiPdf = ROOT.RooProdPdf("sigD0_sigJpsiPdf", "Jpsi Pdf * D0 Pdf", ROOT.RooArgList(cbPdfD0, cbPdfJpsi))
    sigD0_sigPsi2SPdf = ROOT.RooProdPdf("sigD0_sigPsi2SPdf", "Psi2S Pdf * D0 Pdf", ROOT.RooArgList(cbPdfD0, cbPdfPsi2S))
    bkgD0_sigJpsiPdf = ROOT.RooProdPdf("bkgD0_sigJpsiPdf", "Jpsi Pdf * Bkg2 Pdf", ROOT.RooArgList(chebyPdfD0, cbPdfJpsi))
    bkgD0_sigPsi2SPdf = ROOT.RooProdPdf("bkgD0_sigPsi2SPdf", "Psi2S Pdf * Bkg2 Pdf", ROOT.RooArgList(chebyPdfD0, cbPdfPsi2S))
    reflD0_sigJpsiPdf = ROOT.RooProdPdf("reflD0_sigJpsiPdf", "Jpsi Pdf * refl Pdf", ROOT.RooArgList(templatePdfReflD0, cbPdfJpsi))
    reflD0_sigPsi2SPdf = ROOT.RooProdPdf("reflD0_sigPsi2SPdf", "Psi2S Pdf * refl Pdf", ROOT.RooArgList(templatePdfReflD0, cbPdfPsi2S))
    
    bkgJpsi_sigD0Pdf = ROOT.RooProdPdf("bkgJpsi_sigD0Pdf", "D0 Pdf * Bkg1 Pdf", ROOT.RooArgList(chebyPdfJpsi, cbPdfD0))
    bkgD0_bkgJpsiPdf = ROOT.RooProdPdf("bkgD0_bkgJpsiPdf", "Bkg1 Pdf * Bkg2 Pdf", ROOT.RooArgList(chebyPdfD0, chebyPdfJpsi))
    reflD0_bkgJpsiPdf = ROOT.RooProdPdf("reflD0_bkgJpsiPdf", "Bkg1 Pdf * refl Pdf", ROOT.RooArgList(templatePdfReflD0, chebyPdfJpsi))

    # Sum all Pdfs
    model = ROOT.RooAddPdf("model", "sigD0_sigJpsiPdf + bkgD0_sigJpsiPdf + reflD0_sigJpsiPdf + sigD0_sigPsi2SPdf + bkgD0_sigPsi2SPdf + reflD0_sigPsi2SPdf + bkgJpsi_sigD0Pdf + bkgD0_bkgJpsiPdf + reflD0_bkgJpsiPdf", ROOT.RooArgList(sigD0_sigJpsiPdf, bkgD0_sigJpsiPdf, reflD0_sigJpsiPdf, sigD0_sigPsi2SPdf, bkgD0_sigPsi2SPdf, reflD0_sigPsi2SPdf,bkgJpsi_sigD0Pdf, bkgD0_bkgJpsiPdf, reflD0_bkgJpsiPdf), ROOT.RooArgList(nJPsiD0, nBkgJPsi, nReflJPsi, nPsi2SD0, nBkgPsi2S, nReflPsi2S,nBkgD0, nBkgBkg, nReflBkg))

    if config["fit"]["toy_mc"]:
        # Generate toy sample
        data  = sigD0_sigJpsiPdf.generate(ROOT.RooArgSet(mD0, mJpsi), genJpsiD0)
        sig1bkg1Sample  = bkgD0_sigJpsiPdf.generate(ROOT.RooArgSet(mD0, mJpsi), genBkgJpsi)
        sig2bkg2Sample  = bkgJpsi_sigD0Pdf.generate(ROOT.RooArgSet(mJpsi, mD0), genBkgD0)
        bkg1bkg2Sample = bkgD0_bkgJpsiPdf.generate(ROOT.RooArgSet(mD0, mJpsi), genBkgBkg)

        data.append(sig1bkg1Sample)
        data.append(sig2bkg2Sample)
        data.append(bkg1bkg2Sample)

        if not config["fit"]["unbinned"]:
            sample = data.createHistogram("data m_{J/#psi},m_{D0}", mD0, Binning=50, YVar=dict(var=mJpsi, Binning=50))
    else:
        if config["fit"]["unbinned"]:
            fIn = ROOT.TFile(config["inputs"]["data"], "READ")
            sample = fIn.Get(config["inputs"]["tree"])
        else:
            fIn = ROOT.TFile(config["inputs"]["data"], "READ")
            sample = fIn.Get(config["inputs"]["hist"])
    if config["fit"]["unbinned"]:
        sampleToFit = ROOT.RooDataSet("dataTree", "dataTree", [mD0, mJpsi], Import=sample)
    else:
        sampleToFit = ROOT.RooDataHist("dataHist", "dataHist", [mD0, mJpsi], Import=sample)

    # Fit
    fitResult = model.fitTo(sampleToFit, ROOT.RooFit.PrintLevel(3), ROOT.RooFit.Optimize(1), ROOT.RooFit.Hesse(1), ROOT.RooFit.Minos(1), ROOT.RooFit.Strategy(2), ROOT.RooFit.Save(1))

    modelHist = model.createHistogram("model m_{D0}, m_{J/#psi}", mD0, Binning=50, YVar=dict(var=mJpsi, Binning=50))
    modelHist.SetLineColor(ROOT.kRed)
    modelHist.SetLineWidth(1)

    mJpsiframe = mJpsi.frame(Title="Dimuon")
    sampleToFit.plotOn(mJpsiframe, ROOT.RooFit.Binning(50))
    model.plotOn(mJpsiframe)
    model.plotOn(mJpsiframe, Name={"sigD0_sigJpsiPdf"}, Components={sigD0_sigJpsiPdf}, LineStyle="--", LineColor=ROOT.kRed+1)
    model.plotOn(mJpsiframe, Name={"sigD0_sigPsi2SPdf"}, Components={sigD0_sigPsi2SPdf}, LineStyle=5, LineColor=ROOT.kRed+1)
    model.plotOn(mJpsiframe, Name={"bkgD0_sigJpsiPdf"}, Components={bkgD0_sigJpsiPdf}, LineStyle="--", LineColor=ROOT.kAzure+4)
    model.plotOn(mJpsiframe, Name={"bkgD0_sigPsi2SPdf"}, Components={bkgD0_sigPsi2SPdf}, LineStyle=5, LineColor=ROOT.kAzure+4)
    model.plotOn(mJpsiframe, Name={"bkgJpsi_sigD0Pdf"}, Components={bkgJpsi_sigD0Pdf}, LineStyle="--", LineColor=ROOT.kGreen+1)
    model.plotOn(mJpsiframe, Name={"bkgD0_bkgJpsiPdf"}, Components={bkgD0_bkgJpsiPdf}, LineStyle="--", LineColor=ROOT.kMagenta)
    model.plotOn(mJpsiframe, Name={"reflD0_sigJpsiPdf"}, Components={reflD0_sigJpsiPdf}, LineStyle="--", LineColor=ROOT.kYellow+2)
    model.plotOn(mJpsiframe, Name={"reflD0_sigPsi2SPdf"}, Components={reflD0_sigPsi2SPdf}, LineStyle=5, LineColor=ROOT.kYellow+2)
    model.plotOn(mJpsiframe, Name={"reflD0_bkgJpsiPdf"}, Components={reflD0_bkgJpsiPdf}, LineStyle="--", LineColor=ROOT.kGray+1)

    mD0frame = mD0.frame(Title="D-meson")
    sampleToFit.plotOn(mD0frame, ROOT.RooFit.Binning(50))
    model.plotOn(mD0frame)
    model.plotOn(mD0frame, Name={"sigD0_sigJpsiPdf"}, Components={sigD0_sigJpsiPdf}, LineStyle="--", LineColor=ROOT.kRed+1)
    model.plotOn(mD0frame, Name={"sigD0_sigPsi2SPdf"}, Components={sigD0_sigPsi2SPdf}, LineStyle=5, LineColor=ROOT.kRed+1)
    model.plotOn(mD0frame, Name={"bkgD0_sigJpsiPdf"}, Components={bkgD0_sigJpsiPdf}, LineStyle="--", LineColor=ROOT.kAzure+4)
    model.plotOn(mD0frame, Name={"bkgD0_sigPsi2SPdf"}, Components={bkgD0_sigPsi2SPdf}, LineStyle=5, LineColor=ROOT.kAzure+4)
    model.plotOn(mD0frame, Name={"bkgJpsi_sigD0Pdf"}, Components={bkgJpsi_sigD0Pdf}, LineStyle="--", LineColor=ROOT.kGreen+1)
    model.plotOn(mD0frame, Name={"bkgD0_bkgJpsiPdf"}, Components={bkgD0_bkgJpsiPdf}, LineStyle="--", LineColor=ROOT.kMagenta)
    model.plotOn(mD0frame, Name={"reflD0_sigJpsiPdf"}, Components={reflD0_sigJpsiPdf}, LineStyle="--", LineColor=ROOT.kYellow+2)
    model.plotOn(mD0frame, Name={"reflD0_sigPsi2SPdf"}, Components={reflD0_sigPsi2SPdf}, LineStyle=5, LineColor=ROOT.kYellow+2)
    model.plotOn(mD0frame, Name={"reflD0_bkgJpsiPdf"}, Components={reflD0_bkgJpsiPdf}, LineStyle="--", LineColor=ROOT.kGray+1)

    legend = ROOT.TLegend(0.59, 0.59, 0.89, 0.89)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(mJpsiframe.findObject("sigD0_sigJpsiPdf"), "sig. J/#psi - sig. D0", "L")
    legend.AddEntry(mJpsiframe.findObject("sigD0_sigPsi2SPdf"), "sig. #psi(2S) - sig. D0", "L")
    legend.AddEntry(mJpsiframe.findObject("bkgD0_sigJpsiPdf"), "sig. J/#psi - bkg. D0", "L")
    legend.AddEntry(mJpsiframe.findObject("bkgD0_sigPsi2SPdf"), "sig. #psi(2S) - bkg. D0", "L")
    legend.AddEntry(mJpsiframe.findObject("reflD0_sigJpsiPdf"), "sig. J/#psi - refl. D0", "L")
    legend.AddEntry(mJpsiframe.findObject("reflD0_sigPsi2SPdf"), "sig. #psi(2S) - refl. D0", "L")
    legend.AddEntry(mJpsiframe.findObject("bkgJpsi_sigD0Pdf"), "bkg. #psi - sig. D0", "L")
    legend.AddEntry(mJpsiframe.findObject("bkgD0_bkgJpsiPdf"), "bkg. #psi - bkg. D0", "L")
    legend.AddEntry(mJpsiframe.findObject("reflD0_bkgJpsiPdf"), "bkg. #psi - refl. D0", "L")
    legend.Draw()

    canvasFitJpsi = ROOT.TCanvas("canvasFitJpsi", "canvasFitJpsi", 800, 600)
    mJpsiframe.GetYaxis().SetTitleOffset(1.4)
    mJpsiframe.Draw()
    legend.Draw()
    canvasFitJpsi.Update()
    canvasFitJpsi.SaveAs(f'{config["output"]["figures"]}/projected_jpsi_fit.pdf')
    
    canvasFitJpsi.SetLogy()
    mJpsiframe.GetYaxis().SetRangeUser(1, sampleToFit.numEntries())
    canvasFitJpsi.Update()
    canvasFitJpsi.SaveAs(f'{config["output"]["figures"]}/projected_jpsi_fit_logy.pdf')

    canvasFitD0 = ROOT.TCanvas("canvasFitD0", "canvasFitD0", 800, 600)
    mD0frame.GetYaxis().SetTitleOffset(1.4)
    mD0frame.Draw()
    legend.Draw()
    canvasFitD0.Update()
    canvasFitD0.SaveAs(f'{config["output"]["figures"]}/projected_d0_fit.pdf')

    #canvasFit = ROOT.TCanvas("canvasFit", "canvasFit", 1200, 600)
    #canvasFit.Divide(2, 1)

    #canvasFit.cd(1)
    #mJpsiframe.GetYaxis().SetTitleOffset(1.4)
    #mJpsiframe.Draw()
    #legend.Draw()

    #canvasFit.cd(2)
    #mD0frame.GetYaxis().SetTitleOffset(1.4)
    #mD0frame.Draw()
    #legend.Draw()
    #canvasFit.Update()
    #canvasFit.SaveAs(f'{config["output"]["figures"]}/projected_fit.pdf')

    canvasFitHist3D = ROOT.TCanvas("canvasFitHist3D", "canvasFitHist3D", 1000, 1000)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gPad.SetLeftMargin(0.15)
    
    if config["fit"]["unbinned"]:
        sample.Draw("fMass : fMassD0 >> hist_fMass_fMassD0", f'fMass > {minFitRangeJpsi} && fMass < {maxFitRangeJpsi} && fMassD0 > {minFitRangeD0} && fMassD0 < {maxFitRangeD0}', "LEGO2")
        htemp = ROOT.gPad.GetPrimitive("hist_fMass_fMassD0")
        htemp.GetXaxis().SetTitleOffset(2.0)
        htemp.GetYaxis().SetTitleOffset(2.0)
        htemp.GetZaxis().SetTitleOffset(2.0)
        htemp.GetXaxis().SetTitle("#it{m}_{#pi#it{K}} (GeV/#it{c}^{2})")
        htemp.GetYaxis().SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})")
        htemp.Rebin2D(2, 2)
        htemp.Draw("LEGO2")
    else:
        sample.GetXaxis().SetTitleOffset(2.0)
        sample.GetYaxis().SetTitleOffset(2.0)
        sample.GetZaxis().SetTitleOffset(2.0)
        sample.Draw("LEGO2")
    
    modelHist.Draw("SURF SAME")
    canvasFitHist3D.Update()

    fitResult.Print()

    fOut = ROOT.TFile(f'{config["output"]["directory"]}/myTest.root', "RECREATE")
    #canvasFit.Write()
    canvasFitJpsi.Write()
    canvasFitD0.Write()
    canvasFitHist3D.Write()
    fOut.Close()

    input()
    exit()

def upper_limit():
    toy_mc = False
    # Variables
    mD0   = ROOT.RooRealVar("D-meson mass", "#it{m}_{#pi#it{K}} (GeV/#it{c}^{2})", 1.75, 2.00)
    mJpsi = ROOT.RooRealVar("Dimuon mass", "#it{m}_{#mu#mu} (GeV/#it{c}^{2})", 2.50, 4.00)

    # Yields parameters
    genJpsiD0  = 5
    genBkgJpsi = 100
    genBkgD0   = 100
    genBkgBkg  = 500

    if toy_mc:
        nJPsiD0  = ROOT.RooRealVar("nJPsiD0", "number of JPsi-D0", genJpsiD0, 0., 50.)
        nBkgJPsi = ROOT.RooRealVar("nBkgJPsi", "number of JPsi-Bkg", genBkgJpsi, 0., 200.)
        nBkgD0   = ROOT.RooRealVar("nBkgD0", "number of D0-Bkg", genBkgD0, 0., 200.)
        nBkgBkg  = ROOT.RooRealVar("nBkgBkg", "number of Bkg-Bkg", genBkgBkg, 0., 500)

    else:
        nJPsiD0  = ROOT.RooRealVar("nJPsiD0", "number of JPsi-D0", genJpsiD0, 0., genJpsiD0 * 10.)
        nBkgJPsi = ROOT.RooRealVar("nBkgJPsi", "number of JPsi-Bkg", genBkgJpsi, 0., genBkgJpsi * 10.)
        nBkgD0   = ROOT.RooRealVar("nBkgD0", "number of D0-Bkg", genBkgD0, 0., genBkgD0 * 10.)
        nBkgBkg  = ROOT.RooRealVar("nBkgBkg", "number of Bkg-Bkg", genBkgBkg, 0., genBkgBkg * 10.)

    # Pdfs
    meanJpsi  = ROOT.RooRealVar("meanJpsi", "meanJpsi", 3.0741e+00)
    sigmaJpsi = ROOT.RooRealVar("sigmaJpsi", "sigmaJpsi", 9.5902e-02)
    alphaJpsi = ROOT.RooRealVar("alphaJpsi", "alphaJpsi", 1.0632e+00)
    nJpsi     = ROOT.RooRealVar("nJpsi", "nJpsi", 1.0000e+01)
    cbPdfJpsi = ROOT.RooCBShape("cbPdfJpsi", "Crystal Ball J/psi", mJpsi, meanJpsi, sigmaJpsi, alphaJpsi, nJpsi)

    parsJpsi = [-9.8794e-01, 3.0372e-01, -8.5216e-02]
    chebyParsJpsi = [ROOT.RooRealVar(f"cheb_coeff_jpsi_{i}", f"Coeff_jpsi_{i}", parsJpsi[i]) for i in range(3)]
    chebyPdfJpsi = ROOT.RooChebychev("chebyPdfJpsi", "Cheby for Bkg1", mJpsi, ROOT.RooArgList(*chebyParsJpsi))

    meanD0  = ROOT.RooRealVar("meanD0", "meanD0", 1.8600e+00)
    sigmaD0 = ROOT.RooRealVar("sigmaD0", "sigmaD0", 2.2635e-02)
    gausPdfD0 = ROOT.RooGaussian("gausPdfD0", "Gaussian D0", mD0, meanD0, sigmaD0)

    parsD0 = [-1.7735e-01, -2.5525e-02, 1.6542e-03]
    chebyParsD0 = [ROOT.RooRealVar(f"cheb_coeff_d0_{i}", f"Coeff_d0_{i}", parsD0[i]) for i in range(3)]
    chebyPdfD0 = ROOT.RooChebychev("chebyPdfD0", "Cheby for Bkg2", mD0, ROOT.RooArgList(*chebyParsD0))

    # Product of Pdfs
    sigD0_sigJpsiPdf = ROOT.RooProdPdf("sigD0_sigJpsiPdf", "Jpsi Pdf * D0 Pdf", ROOT.RooArgList(gausPdfD0, cbPdfJpsi))
    bkgD0_sigJpsiPdf = ROOT.RooProdPdf("bkgD0_sigJpsiPdf", "Jpsi Pdf * Bkg2 Pdf", ROOT.RooArgList(chebyPdfD0, cbPdfJpsi))
    bkgJpsi_sigD0Pdf = ROOT.RooProdPdf("bkgJpsi_sigD0Pdf", "D0 Pdf * Bkg1 Pdf", ROOT.RooArgList(chebyPdfJpsi, gausPdfD0))
    bkgD0_bkgJpsiPdf = ROOT.RooProdPdf("bkgD0_bkgJpsiPdf", "Bkg1 Pdf * Bkg2 Pdf", ROOT.RooArgList(chebyPdfD0, chebyPdfJpsi))

    # Sum all Pdfs
    model = ROOT.RooAddPdf("model", "sigD0_sigJpsiPdf + bkgD0_sigJpsiPdf + bkgJpsi_sigD0Pdf + bkgD0_bkgJpsiPdf", ROOT.RooArgList(sigD0_sigJpsiPdf, bkgD0_sigJpsiPdf, bkgJpsi_sigD0Pdf, bkgD0_bkgJpsiPdf), ROOT.RooArgList(nJPsiD0, nBkgJPsi, nBkgD0, nBkgBkg))

    if toy_mc:
        # Generate toy sample
        data  = sigD0_sigJpsiPdf.generate(ROOT.RooArgSet(mD0, mJpsi), genJpsiD0)
        sig1bkg1Sample  = bkgD0_sigJpsiPdf.generate(ROOT.RooArgSet(mD0, mJpsi), genBkgJpsi)
        sig2bkg2Sample  = bkgJpsi_sigD0Pdf.generate(ROOT.RooArgSet(mJpsi, mD0), genBkgD0)
        bkg1bkg2Sample = bkgD0_bkgJpsiPdf.generate(ROOT.RooArgSet(mD0, mJpsi), genBkgBkg)

        data.append(sig1bkg1Sample)
        data.append(sig2bkg2Sample)
        data.append(bkg1bkg2Sample)

        histJpsiD0 = data.createHistogram("data m_{J/#psi},m_{D0}", mD0, Binning=50, YVar=dict(var=mJpsi, Binning=50))
    else:
        fIn = ROOT.TFile("data/histograms.root", "READ")
        histJpsiD0= fIn.Get("hSparseJPsiDmeson_proj_3_2")

    dataHist = ROOT.RooDataHist("dataHist", "dataHist", [mD0, mJpsi], Import=histJpsiD0)

    nuisanceParameters = {meanJpsi, sigmaJpsi, alphaJpsi, meanD0, sigmaD0}
    wspace = ROOT.RooWorkspace()
    modelConfig = ROOT.RooStats.ModelConfig(wspace)
    modelConfig.SetPdf(model)
    modelConfig.SetParametersOfInterest({nJPsiD0})
    modelConfig.SetNuisanceParameters(nuisanceParameters)
    #modelConfig.SetGlobalObservables({nBkgBkg, nBkgD0, nBkgJPsi})
    modelConfig.SetObservables({mD0, mJpsi})
    modelConfig.SetName("ModelConfig")
    wspace.Import(modelConfig)
    wspace.Import(dataHist)
    wspace.SetName("w")

    # First, let's use a Calculator based on the Profile Likelihood Ratio
    plc = ROOT.RooStats.ProfileLikelihoodCalculator(dataHist, modelConfig)
    #plc.SetTestSize(0.05)
    plc.SetConfidenceLevel(0.95)
    lrinterval = plc.GetInterval()
    
    # Let's make a plot
    dataCanvas = ROOT.TCanvas("dataCanvas")

    plotInt = ROOT.RooStats.LikelihoodIntervalPlot(lrinterval)
    plotInt.SetTitle("Profile Likelihood Ratio and Posterior for S")
    plotInt.Draw()

    # Second, use a Calculator based on the Feldman Cousins technique
    fc = ROOT.RooStats.FeldmanCousins(dataHist, modelConfig)
    fc.UseAdaptiveSampling(True)
    fc.FluctuateNumDataEntries(True)  # number counting analysis: dataset always has 1 entry with N events observed
    fc.SetNBins(20)  # number of points to test per parameter
    fc.SetTestSize(0.05)
    # fc.SaveBeltToFile(True) # optional
    fcint = fc.GetInterval()

    # Get Lower and Upper limits from Profile Calculator
    print("Profile lower limit on s = ", lrinterval.LowerLimit(nJPsiD0))
    print("Profile upper limit on s = ", lrinterval.UpperLimit(nJPsiD0))

    # Get Lower and Upper limits from FeldmanCousins with profile construction
    if fcint:
        fcul = fcint.UpperLimit(nJPsiD0)
        fcll = fcint.LowerLimit(nJPsiD0)
        print("FC lower limit on s = ", fcll)
        print("FC upper limit on s = ", fcul)
        fcllLine = ROOT.TLine(fcll, 0, fcll, 1)
        fculLine = ROOT.TLine(fcul, 0, fcul, 1)
        fcllLine.SetLineColor(ROOT.kRed)
        fculLine.SetLineColor(ROOT.kRed)
        fcllLine.Draw("same")
        fculLine.Draw("same")
        dataCanvas.Update()

    dataCanvas.SaveAs("figures/upper_limit.pdf")
    input()
    exit()

def plot_results():
    LoadStyle()

    letexTitle = ROOT.TLatex()
    letexTitle.SetTextSize(0.05)
    letexTitle.SetNDC()
    letexTitle.SetTextFont(42)

    fIn = ROOT.TFile("outputs/myTest.root", "READ")
    canvasInJpsi = fIn.Get("canvasFitJpsi")
    canvasInD0 = fIn.Get("canvasFitD0")
    listOfPrimitivesJpsi = canvasInJpsi.GetListOfPrimitives()
    listOfPrimitivesD0 = canvasInD0.GetListOfPrimitives()

    frameJpsi = listOfPrimitivesJpsi.At(1)
    frameJpsi.SetTitle(" ")
    frameJpsi.GetXaxis().SetRangeUser(2.55, 4.00)
    frameJpsi.GetXaxis().SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})")
    frameJpsi.GetXaxis().SetTitleOffset(1.1)
    frameJpsi.GetXaxis().SetTitleSize(0.05)
    frameJpsi.GetXaxis().SetLabelSize(0.045)
    frameJpsi.GetYaxis().SetRangeUser(0, 800)
    frameJpsi.GetYaxis().SetTitle("Counts")
    frameJpsi.GetYaxis().SetTitleOffset(1.3)
    frameJpsi.GetYaxis().SetTitleSize(0.05)
    frameJpsi.GetYaxis().SetLabelSize(0.045)

    # Histograms
    histDataJpsi = listOfPrimitivesJpsi.At(2)

    # PDFs
    pdfJpsi = listOfPrimitivesJpsi.At(3)
    pdfJpsiS1S2 = listOfPrimitivesJpsi.At(4)
    pdfJpsiS1B2 = listOfPrimitivesJpsi.At(5)
    pdfJpsiB1S2 = listOfPrimitivesJpsi.At(6)
    pdfJpsiB1B2 = listOfPrimitivesJpsi.At(7)

    pdfJpsiS1S2.SetLineColor(ROOT.kRed+1)
    pdfJpsiS1B2.SetLineColor(ROOT.kAzure+4)
    pdfJpsiB1S2.SetLineColor(ROOT.kGreen+1)
    pdfJpsiB1B2.SetLineColor(ROOT.kOrange+7)

    canvasOutJpsi = TCanvas("canvasOutJpsi", "canvasOutJpsi", 600, 600)
    canvasOutJpsi.SetTickx(1)
    canvasOutJpsi.SetTicky(1)
    frameJpsi.Draw()
    histDataJpsi.Draw("EP SAME")
    pdfJpsi.Draw("SAME")
    pdfJpsiS1S2.Draw("SAME")
    pdfJpsiS1B2.Draw("SAME")
    pdfJpsiB1S2.Draw("SAME")
    pdfJpsiB1B2.Draw("SAME")

    legend1 = TLegend(0.20, 0.72, 0.55, 0.90, " ", "brNDC")
    SetLegend(legend1)
    legend1.SetTextSize(0.045)
    legend1.AddEntry(histDataJpsi, "Data", "EP")
    legend1.AddEntry(pdfJpsi, "Total fit", "L")
    legend1.Draw()

    legend2 = TLegend(0.64, 0.62, 0.89, 0.90, " ", "brNDC")
    SetLegend(legend2)
    legend2.SetTextSize(0.045)
    legend2.AddEntry(pdfJpsiS1S2, "J/#psi + D^{0}", "L")
    legend2.AddEntry(pdfJpsiS1B2, "J/#psi + bkg", "L")
    legend2.AddEntry(pdfJpsiB1S2, "Bkg + D^{0}", "L")
    legend2.AddEntry(pdfJpsiB1B2, "Bkg + bkg", "L")
    legend2.Draw()

    letexTitle.DrawLatex(0.20, 0.88, "pp, #sqrt{#it{s}} = 13.6 TeV, 6#times 10^{11} events")

    canvasOutJpsi.Update()
    canvasOutJpsi.SaveAs("figures/fit_jpsi_projection.pdf")

    frameD0 = listOfPrimitivesD0.At(1)
    frameD0.SetTitle(" ")
    frameD0.GetXaxis().SetRangeUser(1.70, 2.05)
    frameD0.GetXaxis().SetTitle("#it{m}_{#pi#it{K}} (GeV/#it{c}^{2})")
    frameD0.GetXaxis().SetTitleOffset(1.1)
    frameD0.GetXaxis().SetTitleSize(0.05)
    frameD0.GetXaxis().SetLabelSize(0.045)
    frameD0.GetYaxis().SetRangeUser(0, 800)
    frameD0.GetYaxis().SetTitle("Counts")
    frameD0.GetYaxis().SetTitleOffset(1.3)
    frameD0.GetYaxis().SetTitleSize(0.05)
    frameD0.GetYaxis().SetLabelSize(0.045)

    # Histograms
    histDataD0 = listOfPrimitivesD0.At(2)

    # PDFs
    pdfD0 = listOfPrimitivesD0.At(3)
    pdfD0S1S2 = listOfPrimitivesD0.At(4)
    pdfD0S1B2 = listOfPrimitivesD0.At(5)
    pdfD0B1S2 = listOfPrimitivesD0.At(6)
    pdfD0B1B2 = listOfPrimitivesD0.At(7)

    pdfD0S1S2.SetLineColor(ROOT.kRed+1)
    pdfD0S1B2.SetLineColor(ROOT.kAzure+4)
    pdfD0B1S2.SetLineColor(ROOT.kGreen+1)
    pdfD0B1B2.SetLineColor(ROOT.kOrange+7)

    canvasOutD0 = TCanvas("canvasOutD0", "canvasOutD0", 600, 600)
    canvasOutD0.SetTickx(1)
    canvasOutD0.SetTicky(1)
    frameD0.Draw()
    histDataD0.Draw("EP SAME")
    pdfD0.Draw("SAME")
    pdfD0S1S2.Draw("SAME")
    pdfD0S1B2.Draw("SAME")
    pdfD0B1S2.Draw("SAME")
    pdfD0B1B2.Draw("SAME")

    legend1 = TLegend(0.20, 0.72, 0.55, 0.90, " ", "brNDC")
    SetLegend(legend1)
    legend1.SetTextSize(0.045)
    legend1.AddEntry(histDataD0, "Data", "EP")
    legend1.AddEntry(pdfD0, "Total fit", "L")
    legend1.Draw()

    legend2 = TLegend(0.64, 0.62, 0.89, 0.90, " ", "brNDC")
    SetLegend(legend2)
    legend2.SetTextSize(0.045)
    legend2.AddEntry(pdfD0S1S2, "J/#psi + D^{0}", "L")
    legend2.AddEntry(pdfD0S1B2, "J/#psi + bkg", "L")
    legend2.AddEntry(pdfD0B1S2, "Bkg + D^{0}", "L")
    legend2.AddEntry(pdfD0B1B2, "Bkg + bkg", "L")
    legend2.Draw()

    letexTitle.DrawLatex(0.20, 0.88, "pp, #sqrt{#it{s}} = 13.6 TeV, 6#times 10^{11} events")

    canvasOutD0.Update()
    canvasOutD0.SaveAs("figures/fit_d0_projection.pdf")

    input()

if __name__ == '__main__':
    main()
