import sys
import argparse
import yaml
import ROOT
import numpy as np
sys.path.append('../utils')
from plot_library import LoadStyle, SetGraStat, SetGraSyst, SetLegend

def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='config.yml', help='config file name')
    parser.add_argument("--do_prefit", help="Do the pre-fit to data", action="store_true")
    parser.add_argument("--do_weightdata", help="Weight unbinned data sample with Axe", action="store_true")
    parser.add_argument("--do_fit", help="Do the fit to data / toy MC", action="store_true")
    parser.add_argument("--do_upper_limit", help="Do the calculation of the upper limit", action="store_true")
    parser.add_argument("--do_prefilter", help="Apply selections", action="store_true")
    parser.add_argument("--do_plot_results", help="Plot results", action="store_true")
    args = parser.parse_args()
    

    with open(args.cfgFileName, 'r') as yml_cfg:
        inputCfg = yaml.load(yml_cfg, yaml.FullLoader)

    if args.do_prefit:
        prefit()

    if args.do_weightdata:
        weightdata(inputCfg)

    if args.do_fit:
        fit(inputCfg)

    if args.do_upper_limit:
        upper_limit(inputCfg)

    if args.do_prefilter:
        prefilter(inputCfg)

    if args.do_plot_results:
        plot_results(inputCfg)

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

def weightdata(config):
    """
    function to assing weights to the data (see config_fit.yml)
    """

    fIn = ROOT.TFile(config["inputs"]["data"], "READ")
    treeIn = fIn.Get(config["inputs"]["tree"])

    fOutName = config["inputs"]["data"]
    fOut = ROOT.TFile(fOutName.replace(".root", "_weighted.root"), "RECREATE")
    treeOut = treeIn.CloneTree(0)

    fInAxeJpsi = ROOT.TFile(config["inputs"]["axeJpsi"], "READ")
    histAxeJpsi = fInAxeJpsi.Get("histAxeJpsi")

    fInAxeD0 = ROOT.TFile(config["inputs"]["axeD0"], "READ")
    histAxeD0 = fInAxeD0.Get("histAxeD0_prompt")

    weight = np.zeros(1, dtype=np.double)
    treeOut.Branch("weight", weight, "weight/D")

    for entry in range(treeIn.GetEntries()):
        treeIn.GetEntry(entry)
        
        if (abs(treeIn.fRapDmes) > 0.6 or treeIn.fPtDmes > 30): continue
        if (abs(treeIn.fRapJpsi) > 4 or abs(treeIn.fRapJpsi) < 2.5 or treeIn.fPtJpsi > 20): continue

        ptBinD0 = histAxeD0.GetXaxis().FindBin(treeIn.fPtDmes)
        rapBinD0 = histAxeD0.GetYaxis().FindBin(treeIn.fRapDmes)
        #phiBinD0 = histAxeD0.GetZaxis().FindBin(treeIn.fPhiD0) #WARNING: for the moment (pT, y) correction is applied

        ptBinJpsi = histAxeJpsi.GetXaxis().FindBin(treeIn.fPtJpsi)
        rapBinJpsi = histAxeJpsi.GetYaxis().FindBin(abs(treeIn.fRapJpsi)) #WARNING: in the tree is negative, in Axe is positive
        #phiBinJpsi = histAxeJpsi.GetZaxis().FindBin(treeIn.fPhiJpsi) #WARNING: for the moment (pT, y) correction is applied

        #weightD0 = histAxeD0.GetBinContent(ptBinD0, rapBinD0, phiBinD0) #WARNING: 3D histograms
        #weightJpsi = histAxeJpsi.GetBinContent(ptBinJpsi, rapBinJpsi, phiBinJpsi) #WARNING: 3D histograms
        weightD0 = histAxeD0.GetBinContent(ptBinD0, rapBinD0)
        weightJpsi = histAxeJpsi.GetBinContent(ptBinJpsi, rapBinJpsi)
        print(weightD0, weightJpsi)
        weight[0] = 1. / (weightD0 * weightJpsi)

        treeOut.Fill()

    fOut.cd()
    treeOut.Write()
    fOut.Close()
    fIn.Close()

def prefit():
    """
    function to fit the integrated spectrum for D0 and J/psi to extract shape parameters
    """
    LoadStyle()
    ROOT.gStyle.SetPalette(ROOT.kBird)

    latexTitle = ROOT.TLatex()
    latexTitle.SetTextSize(0.040)
    latexTitle.SetNDC()
    latexTitle.SetTextFont(42)

    # Variables
    mD0   = ROOT.RooRealVar("D-meson mass", "#it{m}_{#piK} (GeV/#it{c}^{2})", 1.76, 2.00)
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
    modelD0.plotOn(mD0frame, Name={"Bkg"}, Components={chebyPdfD0}, LineStyle="--", LineColor=ROOT.kAzure+1)

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
    modelJpsi.plotOn(mJpsiframe, Name={"Bkg"}, Components={chebyPdfJpsi}, LineStyle="--", LineColor=ROOT.kAzure+1)

    canvasFitD0 = ROOT.TCanvas("canvasFitD0", "canvasFitD0", 800, 600)
    canvasFitD0.SetTickx(1)
    canvasFitD0.SetTicky(1)
    mD0frame.GetYaxis().SetTitleOffset(1.4)
    mD0frame.Draw()

    latexTitle.DrawLatex(0.65, 0.74, "#it{N}_{D0} = %1.0f #pm %1.0f" % (nSigD0.getVal(), nSigD0.getError()))
    latexTitle.DrawLatex(0.65, 0.68, "#it{#mu}_{D0} = %4.3f #pm %4.3f" % (meanD0.getVal(), meanD0.getError()))
    latexTitle.DrawLatex(0.65, 0.62, "#it{#sigma}_{D0} = %4.3f #pm %4.3f" % (sigmaD0.getVal(), sigmaD0.getError()))
    
    canvasFitD0.Update()

    canvasFitJpsi = ROOT.TCanvas("canvasFitJpsi", "canvasFitJpsi", 800, 600)
    canvasFitJpsi.SetTickx(1)
    canvasFitJpsi.SetTicky(1)
    mJpsiframe.GetYaxis().SetTitleOffset(1.4)
    mJpsiframe.Draw()

    latexTitle.DrawLatex(0.65, 0.74, "#it{N}_{J/#psi} = %1.0f #pm %1.0f" % (nSigJpsi.getVal(), nSigJpsi.getError()))
    latexTitle.DrawLatex(0.65, 0.68, "#it{#mu}_{J/#psi} = %4.3f #pm %4.3f" % (meanJpsi.getVal(), meanJpsi.getError()))
    latexTitle.DrawLatex(0.65, 0.62, "#it{#sigma}_{J/#psi} = %4.3f #pm %4.3f" % (sigmaJpsi.getVal(), sigmaJpsi.getError()))
    
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

    if config["fit"]["JpsiChannel"] == "Jpsi2ee":
        titleSuffix = "e^{+}e^{-}"
    elif config["fit"]["JpsiChannel"] == "Jpsi2mumu":
        titleSuffix = "#mu^{+}#mu^{-}"
    else:
        print("Error: JpsiChannel not defined in the configuration file.")
        sys.exit(1)
        
    mD0   = ROOT.RooRealVar("fMassDmes", "#it{m}_{#piK} (GeV/#it{c}^{2})", minFitRangeD0, maxFitRangeD0)
    mJpsi = ROOT.RooRealVar("fMass", f"#it{{m}}_{{{titleSuffix}}} (GeV/#it{{c}}^{{2}})", minFitRangeJpsi, maxFitRangeJpsi)
    ptD0   = ROOT.RooRealVar("fPtDmes", "#it{p}_{T,#piK} (GeV/#it{c}^{2})", 0, 100)
    ptJpsi = ROOT.RooRealVar("fPtJpsi", f"#it{{p}}_{{T,{titleSuffix}}} (GeV/#it{{c}}^{{2}})", 0, 100)
    dRap = ROOT.RooRealVar("fDeltaY", f"y_{{{titleSuffix}}} - y_{{#piK}}", -5, -1)
    rapJpsi = ROOT.RooRealVar("fRapJpsi", f"y_{{{titleSuffix}}}", -5, -1)
    if config["fit"]["weighted"]:
        weight = ROOT.RooRealVar("weight", "weight", 0, 1e10)
    
    # Yields parameters
    genJpsiD0  = config["fit"]["norm_par_sig_val"][0]
    genBkgJpsi = config["fit"]["norm_par_sig_val"][1]
    genBkgD0   = config["fit"]["norm_par_sig_val"][2]
    genBkgBkg  = config["fit"]["norm_par_sig_val"][3]
    

    nJPsiD0  = ROOT.RooRealVar("nJPsiD0", "number of JPsi-D0", genJpsiD0, config["fit"]["norm_par_sig_lw_lim"][0], config["fit"]["norm_par_sig_up_lim"][0])
    nBkgJPsi = ROOT.RooRealVar("nBkgJPsi", "number of JPsi-Bkg", genBkgJpsi, config["fit"]["norm_par_sig_lw_lim"][1], config["fit"]["norm_par_sig_up_lim"][1])
    nBkgD0   = ROOT.RooRealVar("nBkgD0", "number of D0-Bkg", genBkgD0, config["fit"]["norm_par_sig_lw_lim"][2], config["fit"]["norm_par_sig_up_lim"][2])
    nBkgBkg  = ROOT.RooRealVar("nBkgBkg", "number of Bkg-Bkg", genBkgBkg, config["fit"]["norm_par_sig_lw_lim"][3], config["fit"]["norm_par_sig_up_lim"][3])
    if config["fit"]["add_psi2s"]:
        rPsi2SJPsi  = ROOT.RooRealVar("rPsi2SJPsi", "Psi2S/JPsi ratio", config["fit"]["norm_par_sig_val"][4], config["fit"]["norm_par_sig_lw_lim"][4], config["fit"]["norm_par_sig_up_lim"][4])
        rPsi2SJPsi.setConstant(config["fit"]["norm_par_sig_is_const"][4])
        nPsi2SD0 = ROOT.RooFormulaVar("nPsi2SD0", "rPsi2SJPsi * nJPsiD0", ROOT.RooArgList(rPsi2SJPsi,nJPsiD0))
        nBkgPsi2S = ROOT.RooFormulaVar("nBkgPsi2S", "rPsi2SJPsi * nBkgJPsi", ROOT.RooArgList(rPsi2SJPsi,nBkgJPsi))
    
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
            sample = data.createHistogram("data m_{J/#psi},m_{D0}", mD0, Binning=config["plot_results"]["dataBins"], YVar=dict(var=mJpsi, Binning=config["plot_results"]["dataBins"]))
    else:
        if config["fit"]["unbinned"]:
            fIn = ROOT.TFile(config["inputs"]["data"], "READ")
            sample = fIn.Get(config["inputs"]["tree"])
        else:
            fIn = ROOT.TFile(config["inputs"]["data"], "READ")
            sample = fIn.Get(config["inputs"]["hist"])
    if config["fit"]["unbinned"]:
        if config["fit"]["weighted"]:
            #sampleToFit = ROOT.RooDataSet("dataTree", "dataTree", [mD0, mJpsi, ptD0, ptJpsi, dRap, rapJpsi, weight], Import=sample, "", "weight")
            sampleToFit = ROOT.RooDataSet("dataTree", "dataset with weights", sample, ROOT.RooArgSet(mD0, mJpsi, ptD0, ptJpsi, dRap, rapJpsi, weight), "", "weight")
        else:
            sampleToFit = ROOT.RooDataSet("dataTree", "dataTree", [mD0, mJpsi, ptD0, ptJpsi, dRap, rapJpsi], Import=sample)
    else:
        sampleToFit = ROOT.RooDataHist("dataHist", "dataHist", [mD0, mJpsi, ptD0, ptJpsi, dRap, rapJpsi], Import=sample)

    if config["fit"]["unbinned"]:
        sampleToFit = sampleToFit.reduce(f"fRapJpsi>{config['fit']['min_jpsi_rap']} && fRapJpsi<{config['fit']['max_jpsi_rap']}")

    # Fit
    fitResult = model.fitTo(sampleToFit, ROOT.RooFit.PrintLevel(3), ROOT.RooFit.Optimize(1), ROOT.RooFit.Hesse(1), ROOT.RooFit.Minos(1), ROOT.RooFit.Strategy(2), ROOT.RooFit.Save(1))
    
    # Create sPlot
    sampleToFit_sPlot = sampleToFit.Clone("sampleToFit_sPlot") # cloning the original dataset
    sData = ROOT.RooStats.SPlot("sData", "An SPlot", sampleToFit_sPlot, model, ROOT.RooArgList(nJPsiD0, nBkgJPsi, nBkgD0, nBkgBkg)) # The line creates an SPlot 'sData' and adds columns to 'sampleToFit_sPlot' that represent the weights for various components of the model
    
    nJPsiD0_sw = ROOT.RooRealVar("nJPsiD0_sw","sWeights for JpsiD0", 0, 1000)
    nBkgJPsi_sw = ROOT.RooRealVar("nBkgJPsi_sw","sWeights for BkgJpsi",0 ,1000)
    nBkgD0_sw = ROOT.RooRealVar("nBkgD0_sw","sWeights for BkgD0", 0 ,1000)
    nBkgBkg_sw = ROOT.RooRealVar("nBkgBkg_sw","sWeights for BkgBkg", 0 ,1000)
    
    dataw_JPsiD0 = ROOT.RooDataSet("dataw_JPsiD0", "dataw_JPsiD0", sampleToFit_sPlot, ROOT.RooArgSet(ptD0, ptJpsi, dRap, nJPsiD0_sw), "", "nJPsiD0_sw")
    dataw_BkgJPsi = ROOT.RooDataSet("dataw_BkgJPsi", "dataw_BkgJPsi", sampleToFit_sPlot, ROOT.RooArgSet(ptD0, ptJpsi, dRap, nBkgJPsi_sw), "", "nBkgJPsi_sw")
    dataw_BkgD0 = ROOT.RooDataSet("dataw_BkgD0", "dataw_BkgD0", sampleToFit_sPlot, ROOT.RooArgSet(ptD0, ptJpsi, dRap, nBkgD0_sw), "", "nBkgD0_sw")
    dataw_BkgBkg = ROOT.RooDataSet("dataw_BkgBkg", "dataw_BkgBkg", sampleToFit_sPlot, ROOT.RooArgSet(ptD0, ptJpsi, dRap, nBkgBkg_sw), "", "nBkgBkg_sw")

    fOut = ROOT.TFile(f'{config["output"]["directory"]}/datasets_splot.root', "RECREATE")
    sampleToFit_sPlot.Write()
    dataw_JPsiD0.Write()
    dataw_BkgJPsi.Write()
    dataw_BkgD0.Write()
    dataw_BkgBkg.Write()
    
    fOut.Close()
    
    #drawing
    modelHist = model.createHistogram("model m_{D0}, m_{J/#psi}", mD0, Binning=config["plot_results"]["dataBins"], YVar=dict(var=mJpsi, Binning=config["plot_results"]["dataBins"]))
    modelHist.SetLineColor(ROOT.kRed)
    modelHist.SetLineWidth(1)

    mJpsiframe = mJpsi.frame(Title=" ")
    sampleToFit.plotOn(mJpsiframe, ROOT.RooFit.Name("data"), ROOT.RooFit.Binning(config["plot_results"]["dataBins"]))
    model.plotOn(mJpsiframe, Name={"model"})
    model.plotOn(mJpsiframe, Name={"sigD0_sigJpsiPdf"}, Components={sigD0_sigJpsiPdf}, LineStyle=5, LineColor=ROOT.kRed+1)
    model.plotOn(mJpsiframe, Name={"bkgD0_sigJpsiPdf"}, Components={bkgD0_sigJpsiPdf}, LineStyle=5, LineColor=ROOT.kAzure+1)
    model.plotOn(mJpsiframe, Name={"bkgJpsi_sigD0Pdf"}, Components={bkgJpsi_sigD0Pdf}, LineStyle=2, LineColor=ROOT.kGreen+1)
    model.plotOn(mJpsiframe, Name={"bkgD0_bkgJpsiPdf"}, Components={bkgD0_bkgJpsiPdf}, LineStyle=2, LineColor=ROOT.kMagenta+1)
    model.plotOn(mJpsiframe, Name={"reflD0_sigJpsiPdf"}, Components={reflD0_sigJpsiPdf}, LineStyle=5, LineColor=ROOT.kYellow+2)
    model.plotOn(mJpsiframe, Name={"reflD0_bkgJpsiPdf"}, Components={reflD0_bkgJpsiPdf}, LineStyle=2, LineColor=ROOT.kGray+1)
    if config["fit"]["add_psi2s"]:
        model.plotOn(mJpsiframe, Name={"sigD0_sigPsi2SPdf"}, Components={sigD0_sigPsi2SPdf}, LineStyle=3, LineColor=ROOT.kRed+1)
        model.plotOn(mJpsiframe, Name={"bkgD0_sigPsi2SPdf"}, Components={bkgD0_sigPsi2SPdf}, LineStyle=3, LineColor=ROOT.kAzure+1)
        model.plotOn(mJpsiframe, Name={"reflD0_sigPsi2SPdf"}, Components={reflD0_sigPsi2SPdf}, LineStyle=3, LineColor=ROOT.kYellow+2)

    mD0frame = mD0.frame(Title=" ")
    sampleToFit.plotOn(mD0frame, ROOT.RooFit.Name("data"), ROOT.RooFit.Binning(config["plot_results"]["dataBins"]))
    model.plotOn(mD0frame, Name={"model"})
    model.plotOn(mD0frame, Name={"sigD0_sigJpsiPdf"}, Components={sigD0_sigJpsiPdf}, LineStyle=5, LineColor=ROOT.kRed+1)
    model.plotOn(mD0frame, Name={"bkgD0_sigJpsiPdf"}, Components={bkgD0_sigJpsiPdf}, LineStyle=5, LineColor=ROOT.kAzure+1)
    model.plotOn(mD0frame, Name={"bkgJpsi_sigD0Pdf"}, Components={bkgJpsi_sigD0Pdf}, LineStyle=2, LineColor=ROOT.kGreen+1)
    model.plotOn(mD0frame, Name={"bkgD0_bkgJpsiPdf"}, Components={bkgD0_bkgJpsiPdf}, LineStyle=2, LineColor=ROOT.kMagenta+1)
    model.plotOn(mD0frame, Name={"reflD0_sigJpsiPdf"}, Components={reflD0_sigJpsiPdf}, LineStyle=5, LineColor=ROOT.kYellow+2)
    model.plotOn(mD0frame, Name={"reflD0_bkgJpsiPdf"}, Components={reflD0_bkgJpsiPdf}, LineStyle=2, LineColor=ROOT.kGray+1)
    if config["fit"]["add_psi2s"]:
        model.plotOn(mD0frame, Name={"sigD0_sigPsi2SPdf"}, Components={sigD0_sigPsi2SPdf}, LineStyle=3, LineColor=ROOT.kRed+1)
        model.plotOn(mD0frame, Name={"bkgD0_sigPsi2SPdf"}, Components={bkgD0_sigPsi2SPdf}, LineStyle=3, LineColor=ROOT.kAzure+1)
        model.plotOn(mD0frame, Name={"reflD0_sigPsi2SPdf"}, Components={reflD0_sigPsi2SPdf}, LineStyle=3, LineColor=ROOT.kYellow+2)
    
    legend_comp = ROOT.TLegend(0.59, 0.57, 0.89, 0.87)
    legend_comp.SetBorderSize(0)
    legend_comp.SetFillStyle(0)
    legend_comp.AddEntry(mD0frame.findObject("data"), "data", "PE")
    legend_comp.AddEntry(mD0frame.findObject("model"), "total fit", "L")

    if config["fit"]["add_psi2s"]:
        hist_Jpsi = ROOT.TH1F("hist_Jpsi","", 100, minFitRangeD0, maxFitRangeD0); hist_Jpsi.SetLineWidth(3); hist_Jpsi.SetLineStyle(5); hist_Jpsi.SetLineColor(ROOT.kBlack)
        hist_psi2s = ROOT.TH1F("hist_psi2s","", 100, minFitRangeD0, maxFitRangeD0); hist_psi2s.SetLineWidth(3); hist_psi2s.SetLineStyle(3); hist_psi2s.SetLineColor(ROOT.kBlack)
        hist_sigD0_sigJpsi = ROOT.TH1F("hist_sigD0_sigJpsi","", 100, minFitRangeD0, maxFitRangeD0); hist_sigD0_sigJpsi.SetLineWidth(3); hist_sigD0_sigJpsi.SetLineStyle(5); hist_sigD0_sigJpsi.SetLineColor(ROOT.kRed+1)
        hist_bkgD0_sigJpsi = ROOT.TH1F("hist_bkgD0_sigJpsi","", 100, minFitRangeD0, maxFitRangeD0); hist_bkgD0_sigJpsi.SetLineWidth(3); hist_bkgD0_sigJpsi.SetLineStyle(5); hist_bkgD0_sigJpsi.SetLineColor(ROOT.kAzure+1)
        hist_bkgJpsi_sigD0 = ROOT.TH1F("hist_bkgJpsi_sigD0","", 100, minFitRangeD0, maxFitRangeD0); hist_bkgJpsi_sigD0.SetLineWidth(3); hist_bkgJpsi_sigD0.SetLineStyle(2); hist_bkgJpsi_sigD0.SetLineColor(ROOT.kGreen+1)
        hist_bkgD0_bkgJpsi = ROOT.TH1F("hist_bkgD0_bkgJpsi","", 100, minFitRangeD0, maxFitRangeD0); hist_bkgD0_bkgJpsi.SetLineWidth(3); hist_bkgD0_bkgJpsi.SetLineStyle(2); hist_bkgD0_bkgJpsi.SetLineColor(ROOT.kMagenta+1)
        hist_reflD0_sigJpsi = ROOT.TH1F("hist_reflD0_sigJpsi","", 100, minFitRangeD0, maxFitRangeD0); hist_reflD0_sigJpsi.SetLineWidth(3); hist_reflD0_sigJpsi.SetLineStyle(5); hist_reflD0_sigJpsi.SetLineColor(ROOT.kYellow+2)
        hist_reflD0_bkgJpsi = ROOT.TH1F("hist_reflD0_bkgJpsi","", 100, minFitRangeD0, maxFitRangeD0); hist_reflD0_bkgJpsi.SetLineWidth(3); hist_reflD0_bkgJpsi.SetLineStyle(2); hist_reflD0_bkgJpsi.SetLineColor(ROOT.kGray+1)
        legend_comp.AddEntry(hist_sigD0_sigJpsi, "sig. #psi - sig. D^{0}", "L")
        legend_comp.AddEntry(hist_bkgD0_sigJpsi, "sig. #psi - bkg. D^{0}", "L")
        legend_comp.AddEntry(hist_reflD0_sigJpsi, "sig. #psi - refl. D^{0}", "L")
        legend_comp.AddEntry(hist_Jpsi,"","")
        legend_comp.AddEntry(hist_psi2s,"","")
        legend_comp.AddEntry(hist_bkgJpsi_sigD0, "bkg. #psi - sig. D^{0}", "L")
        legend_comp.AddEntry(hist_bkgD0_bkgJpsi, "bkg. #psi - bkg. D^{0}", "L")
        legend_comp.AddEntry(hist_reflD0_bkgJpsi, "bkg. #psi - refl. D^{0}", "L")

        legend_comp2 = ROOT.TLegend(legend_comp.GetX1()+0.05, legend_comp.GetY1(), legend_comp.GetX2()+0.05, legend_comp.GetY2())
        legend_comp2.SetBorderSize(0)
        legend_comp2.SetFillStyle(0)
        legend_comp2.AddEntry(hist_sigD0_sigJpsi, "", "") #dummy entries to get the right position
        legend_comp2.AddEntry(hist_sigD0_sigJpsi, "", "")
        legend_comp2.AddEntry(hist_sigD0_sigJpsi, "", "")
        legend_comp2.AddEntry(hist_bkgD0_sigJpsi, "", "")
        legend_comp2.AddEntry(hist_reflD0_sigJpsi, "", "")
        legend_comp2.AddEntry(hist_Jpsi,"J/#psi","L")
        legend_comp2.AddEntry(hist_psi2s,"#psi(2S)","L")
        legend_comp2.AddEntry(hist_bkgJpsi_sigD0, "", "")
        legend_comp2.AddEntry(hist_bkgD0_bkgJpsi, "", "")
        legend_comp2.AddEntry(hist_reflD0_bkgJpsi, "", "")


    else: #if the psi2s is not included just use the same for D0
        legend_comp.AddEntry(mD0frame.findObject("sigD0_sigJpsiPdf"), "sig. J/#psi - sig. D^{0}", "L")
        legend_comp.AddEntry(mD0frame.findObject("bkgD0_sigJpsiPdf"), "sig. J/#psi - bkg. D^{0}", "L")
        legend_comp.AddEntry(mD0frame.findObject("reflD0_sigJpsiPdf"), "sig. J/#psi - refl. D^{0}", "L")
        legend_comp.AddEntry(mD0frame.findObject("bkgJpsi_sigD0Pdf"), "bkg. J/#psi - sig. D^{0}", "L")
        legend_comp.AddEntry(mD0frame.findObject("bkgD0_bkgJpsiPdf"), "bkg. J/#psi - bkg. D^{0}", "L")
        legend_comp.AddEntry(mD0frame.findObject("reflD0_bkgJpsiPdf"), "bkg. J/#psi - refl. D^{0}", "L")

    
    #legend_comp.Draw()

    latexTitle = ROOT.TLatex()
    latexTitle.SetTextSize(0.05)
    latexTitle.SetNDC()
    latexTitle.SetTextFont(42)

    latexRap = ROOT.TLatex()
    latexRap.SetTextSize(0.03)
    latexRap.SetNDC()
    latexRap.SetTextFont(42)

    
    
    canvasFitJpsi = ROOT.TCanvas("canvasFitJpsi", "canvasFitJpsi", 800, 800)
    canvasFitJpsi.SetTickx(1)
    canvasFitJpsi.SetTicky(1)
    mJpsiframe.GetYaxis().SetRangeUser(config["plot_results"]["jpsiFrame"]["y_range"][0], config["plot_results"]["jpsiFrame"]["y_range"][1])
    mJpsiframe.GetYaxis().SetTitleOffset(1.4)
    mJpsiframe.GetYaxis().SetTitle(f'Counts per {round(1000*((config["fit"]["max_fit_range_jpsi"]-config["fit"]["min_fit_range_jpsi"])/config["plot_results"]["dataBins"]))} MeV/#it{{c}}^{{2}}')
    mJpsiframe.Draw()
    legend_comp.Draw()
    if config["fit"]["add_psi2s"]:
        legend_comp2.Draw()
    latexTitle.DrawLatex(0.15, 0.835, "Work in progress")
    latexTitle.DrawLatex(0.15, 0.78, "pp, #sqrt{#it{s}} = 13.6 TeV")
    latexRap.DrawLatex(0.15, 0.72, "|#it{#eta}_{#piK}| < 0.8")
    if config["fit"]["JpsiChannel"] == "Jpsi2mumu":
        latexRap.DrawLatex(0.15, 0.68, f"{config['fit']['min_jpsi_rap']} < #it{{#eta}}_{{#mu#mu}} < {config['fit']['max_jpsi_rap']}") #titleSuffix
    else:
        latexRap.DrawLatex(0.15, 0.68, f"|#it{{#eta}}_{{ee}}| < {config['fit']['max_jpsi_rap']}")
        
    canvasFitJpsi.Update()
    canvasFitJpsi.SaveAs(f'{config["output"]["figures"]}/projected_{config["fit"]["JpsiChannel"]}_jpsi_fit.png')
    
    canvasFitJpsi.SetLogy()
    mJpsiframe.GetYaxis().SetRangeUser(config["plot_results"]["jpsiFrame"]["y_range_log"][0], config["plot_results"]["jpsiFrame"]["y_range_log"][1])
    canvasFitJpsi.Update()
    canvasFitJpsi.SaveAs(f'{config["output"]["figures"]}/projected_{config["fit"]["JpsiChannel"]}_jpsi_fit_logy.png')

    canvasFitD0 = ROOT.TCanvas("canvasFitD0", "canvasFitD0", 800, 800)
    canvasFitD0.SetTickx(1)
    canvasFitD0.SetTicky(1)
    mD0frame.GetYaxis().SetTitleOffset(1.4)
    mD0frame.GetYaxis().SetLabelSize(0.03)
    mD0frame.GetYaxis().SetRangeUser(config["plot_results"]["d0Frame"]["y_range"][0], config["plot_results"]["d0Frame"]["y_range"][1])
    mD0frame.GetYaxis().SetTitle(f'Counts per {round(1000*((config["fit"]["max_fit_range_d0"]-config["fit"]["min_fit_range_d0"])/config["plot_results"]["dataBins"]))} MeV/#it{{c}}^{{2}}')
    mD0frame.Draw()
    legend_comp.Draw()
    if config["fit"]["add_psi2s"]:
        legend_comp2.Draw()
    latexTitle.DrawLatex(0.15, 0.835, "Work in progress")
    latexTitle.DrawLatex(0.15, 0.78, "pp, #sqrt{#it{s}} = 13.6 TeV")
    latexRap.DrawLatex(0.15, 0.72, "|#it{#eta}_{#piK}| < 0.8")
    if config["fit"]["JpsiChannel"] == "Jpsi2mumu":
        latexRap.DrawLatex(0.15, 0.68, f"{config['fit']['min_jpsi_rap']} < #it{{#eta}}_{{#mu#mu}} < {config['fit']['max_jpsi_rap']}") #titleSuffix
    else:
        latexRap.DrawLatex(0.15, 0.68, f"|#it{{#eta}}_{{ee}}| < {config['fit']['max_jpsi_rap']}")

    canvasFitD0.Update()
    canvasFitD0.SaveAs(f'{config["output"]["figures"]}/projected_{config["fit"]["JpsiChannel"]}_d0_fit.png')

    canvasFitD0.SetLogy()
    mD0frame.GetYaxis().SetRangeUser(config["plot_results"]["d0Frame"]["y_range_log"][0], config["plot_results"]["d0Frame"]["y_range_log"][1])
    canvasFitD0.Update()
    canvasFitD0.SaveAs(f'{config["output"]["figures"]}/projected_{config["fit"]["JpsiChannel"]}_d0_fit_logy.png')

    canvasFitHist3D = ROOT.TCanvas("canvasFitHist3D", "canvasFitHist3D", 1000, 1000)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gPad.SetLeftMargin(0.15)
    
    if config["fit"]["unbinned"]:
        sample.Draw("fMass : fMassDmes >> hist_fMass_fMassDmes", f'fMass > {minFitRangeJpsi} && fMass < {maxFitRangeJpsi} && fMassDmes > {minFitRangeD0} && fMassDmes < {maxFitRangeD0}', "LEGO2")
        htemp = ROOT.gPad.GetPrimitive("hist_fMass_fMassDmes")
        htemp.SetTitle(" ")
        htemp.GetXaxis().SetTitleOffset(2.0)
        htemp.GetXaxis().SetLabelSize(0.03)
        htemp.GetYaxis().SetTitleOffset(2.0)
        htemp.GetYaxis().SetLabelSize(0.03)
        htemp.GetZaxis().SetTitleOffset(2.0)
        htemp.GetZaxis().SetLabelSize(0.03)
        htemp.GetXaxis().SetTitle("#it{m}_{#piK} (GeV/#it{c}^{2})")
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
    canvasFitHist3D.SetTitle(" ")
    canvasFitHist3D.SetPhi(230)
    canvasFitHist3D.SetTheta(20)
    
    latexTitle.DrawLatex(0.1, 0.94, "Work in progress")
    latexTitle.DrawLatex(0.1, 0.88, "pp, #sqrt{#it{s}} = 13.6 TeV")
    latexRap.DrawLatex(0.73, 0.9, "|#it{#eta}_{#piK}| < 0.8")
    if config["fit"]["JpsiChannel"] == "Jpsi2mumu":
        latexRap.DrawLatex(0.73, 0.86, f"{config['fit']['min_jpsi_rap']} < #it{{#eta}}_{{#mu#mu}} < {config['fit']['max_jpsi_rap']}") #titleSuffix
    else:
        latexRap.DrawLatex(0.73, 0.86, f"|#it{{#eta}}_{{ee}}| < {config['fit']['max_jpsi_rap']}")

    canvasFitHist3D.Update()
    canvasFitHist3D.SaveAs(f'{config["output"]["figures"]}/fit_2D_{config["fit"]["JpsiChannel"]}_fitSave.png')
    
    fitResult.Print()

    fOut = ROOT.TFile(f'{config["output"]["directory"]}/results.root', "RECREATE")
    #canvasFit.Write()
    canvasFitJpsi.Write()
    canvasFitD0.Write()
    canvasFitHist3D.Write()
    fOut.Close()

    input()
    exit()

def upper_limit(config):
    toy_mc = False
    # Variables
    mD0   = ROOT.RooRealVar("D-meson mass", "#it{m}_{#piK} (GeV/#it{c}^{2})", 1.75, 2.00)
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
        
        histJpsiD0 = data.createHistogram("data m_{J/#psi},m_{D0}", mD0, Binning=config["plot_results"]["dataBins"], YVar=dict(var=mJpsi, Binning=config["plot_results"]["dataBins"]))
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

    dataCanvas.SaveAs(f'{config["output"]["figures"]}/upper_limit.pdf')
    input()
    exit()

def plot_results(config):
    LoadStyle()

    latexTitle = ROOT.TLatex()
    latexTitle.SetTextSize(0.05)
    latexTitle.SetNDC()
    latexTitle.SetTextFont(42)
    
    latexDecay = ROOT.TLatex()
    latexDecay.SetTextSize(0.04)
    latexDecay.SetNDC()
    latexDecay.SetTextFont(42)

    latexRap = ROOT.TLatex()
    latexRap.SetTextSize(0.035)
    latexRap.SetNDC()
    latexRap.SetTextFont(42)

    if config["fit"]["JpsiChannel"] == "Jpsi2mumu":
        decayLabel = "#mu#mu"
    else:
        decayLabel = "ee"
        
    fIn = ROOT.TFile(f'{config["output"]["directory"]}/myTest.root', "READ")
    canvasInJpsi = fIn.Get("canvasFitJpsi")
    canvasInD0 = fIn.Get("canvasFitD0")
    canvasIn3D = fIn.Get("canvasFitHist3D")
    listOfPrimitivesJpsi = canvasInJpsi.GetListOfPrimitives()
    listOfPrimitivesD0 = canvasInD0.GetListOfPrimitives()
    listOfPrimitives3D = canvasIn3D.GetListOfPrimitives()
    
    ## plot Jpsi results
    frameJpsi = listOfPrimitivesJpsi.At(1)
    frameJpsi.SetTitle(" ")
    frameJpsi.GetXaxis().SetRangeUser(config["plot_results"]["jpsiFrame"]["x_range"][0], config["plot_results"]["jpsiFrame"]["x_range"][1])
    frameJpsi.GetXaxis().SetTitle(config["plot_results"]["jpsiFrame"]["x_title"])
    frameJpsi.GetXaxis().SetTitleOffset(config["plot_results"]["jpsiFrame"]["x_title_offset"])
    # frameJpsi.GetXaxis().SetTitleSize(0.05)
    # frameJpsi.GetXaxis().SetLabelSize(0.045)
    frameJpsi.GetYaxis().SetRangeUser(config["plot_results"]["jpsiFrame"]["y_range"][0], config["plot_results"]["jpsiFrame"]["y_range"][1])
    frameJpsi.GetYaxis().SetTitle(f'Counts per {round(1000*((config["fit"]["max_fit_range_jpsi"]-config["fit"]["min_fit_range_jpsi"])/config["plot_results"]["dataBins"]))} MeV/#it{{c}}^{{2}}')
    
    if config["plot_results"]["jpsiFrame"]["y_title"]:
        frameJpsi.GetYaxis().SetTitle(config["plot_results"]["jpsiFrame"]["y_title"])
    frameJpsi.GetYaxis().SetTitleOffset(config["plot_results"]["jpsiFrame"]["y_title_offset"])
    # frameJpsi.GetYaxis().SetTitleSize(0.05)
    # frameJpsi.GetYaxis().SetLabelSize(0.045)

    # Histograms
    histDataJpsi = listOfPrimitivesJpsi.At(2)
    # PDFs
    pdfJpsi = listOfPrimitivesJpsi.At(3)
    pdfJpsiS1S2 = listOfPrimitivesJpsi.At(4)
    pdfJpsiS1B2 = listOfPrimitivesJpsi.At(5)
    pdfJpsiB1S2 = listOfPrimitivesJpsi.At(6)
    pdfJpsiB1B2 = listOfPrimitivesJpsi.At(7)
    pdfJpsiS1F2 = listOfPrimitivesJpsi.At(8)
    pdfJpsiB1F2 = listOfPrimitivesJpsi.At(9)
    if config["fit"]["add_psi2s"]:
        pdfJpsiS2S3 = listOfPrimitivesJpsi.At(10)
        pdfJpsiB2S3 = listOfPrimitivesJpsi.At(11)
        pdfJpsiF2S3 = listOfPrimitivesJpsi.At(12)

    # pdfJpsiS1S2.SetLineColor(ROOT.kRed+1)
    # pdfJpsiS1B2.SetLineColor(ROOT.kAzure+1)
    # pdfJpsiB1S2.SetLineColor(ROOT.kGreen+1)
    # pdfJpsiB1B2.SetLineColor(ROOT.kOrange+7)
    pdfJpsiS1S2.SetLineStyle(5)
    pdfJpsiS1B2.SetLineStyle(5)
    pdfJpsiB1S2.SetLineStyle(2)
    pdfJpsiB1B2.SetLineStyle(2)
    pdfJpsiS1F2.SetLineStyle(5)
    pdfJpsiB1F2.SetLineStyle(2)
    if config["fit"]["add_psi2s"]:
        pdfJpsiS2S3.SetLineStyle(3)
        pdfJpsiB2S3.SetLineStyle(3)
        pdfJpsiF2S3.SetLineStyle(3)

    canvasOutJpsi = ROOT.TCanvas("canvasOutJpsi", "canvasOutJpsi", 800, 800)
    canvasOutJpsi.SetTickx(1)
    canvasOutJpsi.SetTicky(1)
    frameJpsi.Draw("AXIS")
    histDataJpsi.Draw("EP SAME")
    pdfJpsi.Draw("SAME")
    pdfJpsiS1S2.Draw("SAME")
    pdfJpsiS1B2.Draw("SAME")
    pdfJpsiB1S2.Draw("SAME")
    pdfJpsiB1B2.Draw("SAME")
    pdfJpsiS1F2.Draw("SAME")
    pdfJpsiB1F2.Draw("SAME")
    if config["fit"]["add_psi2s"]:
        pdfJpsiS2S3.Draw("SAME")
        pdfJpsiB2S3.Draw("SAME")
        pdfJpsiF2S3.Draw("SAME")

    legend1 = ROOT.TLegend(0.65, 0.65, 0.89, 0.95, " ", "brNDC")
    if config["fit"]["add_psi2s"]:
        legend1 = ROOT.TLegend(0.6, 0.55, 0.8, 0.95, " ", "brNDC")
    SetLegend(legend1)
    legend1.SetTextSize(0.030)
    
    legend1.AddEntry(histDataJpsi, "data", "EP")
    legend1.AddEntry(pdfJpsi, "total fit", "L")

    legend1.AddEntry(pdfJpsiS1S2, "sig. J/#psi - sig. D^{0}", "L")
    legend1.AddEntry(pdfJpsiS1B2, "sig. J/#psi - bkg. #piK", "L")
    legend1.AddEntry(pdfJpsiS1F2, "sig. J/#psi - refl. #piK", "L")
    
    if config["fit"]["add_psi2s"]:
        legend1.AddEntry(pdfJpsiS2S3, "sig. #psi(2S) - sig. D^{0}", "L")
        legend1.AddEntry(pdfJpsiB2S3, "sig. #psi(2S) - bkg. #piK", "L")
        legend1.AddEntry(pdfJpsiF2S3, "sig. #psi(2S) - refl. #piK", "L")
        
        legend1.AddEntry(pdfJpsiB1S2, "bkg. "+decayLabel+" - sig. D^{0}", "L")
        legend1.AddEntry(pdfJpsiB1B2, "bkg. "+decayLabel+" - bkg. #piK", "L")
        legend1.AddEntry(pdfJpsiB1F2, "bkg. "+decayLabel+" - refl. #piK", "L")
    else:
        legend1.AddEntry(pdfJpsiB1S2, "bkg. "+decayLabel+" - sig. D^{0}", "L")
        legend1.AddEntry(pdfJpsiB1B2, "bkg. "+decayLabel+" - bkg. #piK", "L")
        legend1.AddEntry(pdfJpsiB1F2, "bkg. "+decayLabel+" - refl. #piK", "L")
        
    legend1.Draw()

    latexTitle.DrawLatex(0.18, 0.89, "Work in progress")
    latexTitle.DrawLatex(0.18, 0.83, "pp, #sqrt{#it{s}} = 13.6 TeV")
    latexDecay.DrawLatex(0.18, 0.77, "D^{0} #rightarrow #pi^{+}K^{#minus}")
    if config["fit"]["JpsiChannel"] == "Jpsi2mumu":
        latexDecay.DrawLatex(0.18, 0.73, "J/#psi #rightarrow #mu^{+}#mu^{#minus}")
    else:
        latexDecay.DrawLatex(0.18, 0.73, "J/#psi #rightarrow e^{+}e^{#minus}")
    
    latexRap.DrawLatex(0.18, 0.68, "|#it{#eta}_{#piK}| < 0.8")
    if config["fit"]["JpsiChannel"] == "Jpsi2mumu":
        latexRap.DrawLatex(0.18, 0.63, f"{config['fit']['min_jpsi_rap']} < #it{{#eta}}_{{#mu#mu}} < {config['fit']['max_jpsi_rap']}") #titleSuffix
    else:
        latexRap.DrawLatex(0.18, 0.63, f"|#it{{#eta}}_{{ee}}| < {config['fit']['max_jpsi_rap']}")

    canvasOutJpsi.Update()
    canvasOutJpsi.SaveAs(f'{config["output"]["figures"]}/fit_{config["fit"]["JpsiChannel"]}_jpsi_projection.png')
    
    canvasOutJpsi.SetLogy()
    frameJpsi.GetYaxis().SetRangeUser(config["plot_results"]["jpsiFrame"]["y_range_log"][0], config["plot_results"]["jpsiFrame"]["y_range_log"][1])
    frameJpsi.GetYaxis().SetTitleOffset(config["plot_results"]["jpsiFrame"]["y_title_offset_log"])
    canvasOutJpsi.Update()
    canvasOutJpsi.SaveAs(f'{config["output"]["figures"]}/fit_{config["fit"]["JpsiChannel"]}_jpsi_projection_logy.png')

    ## plot D0 results
    frameD0 = listOfPrimitivesD0.At(1)
    frameD0.SetTitle(" ")
    frameD0.GetXaxis().SetRangeUser(config["plot_results"]["d0Frame"]["x_range"][0], config["plot_results"]["d0Frame"]["x_range"][1])
    frameD0.GetXaxis().SetTitle(config["plot_results"]["d0Frame"]["x_title"])
    frameD0.GetXaxis().SetTitleOffset(config["plot_results"]["d0Frame"]["x_title_offset"])
    # frameD0.GetXaxis().SetTitleSize(0.05)
    # frameD0.GetXaxis().SetLabelSize(0.045)
    frameD0.GetYaxis().SetRangeUser(config["plot_results"]["d0Frame"]["y_range"][0], config["plot_results"]["d0Frame"]["y_range"][1])
    frameD0.GetYaxis().SetTitle(f'Counts per {round(1000*((config["fit"]["max_fit_range_d0"]-config["fit"]["min_fit_range_d0"])/config["plot_results"]["dataBins"]))} MeV/#it{{c}}^{{2}}')
    if config["plot_results"]["d0Frame"]["y_title"]:
        frameD0.GetYaxis().SetTitle(config["plot_results"]["d0Frame"]["y_title"])
    frameD0.GetYaxis().SetTitleOffset(config["plot_results"]["d0Frame"]["y_title_offset"])
    # frameD0.GetYaxis().SetTitleSize(0.05)
    # frameD0.GetYaxis().SetLabelSize(0.045)

    # Histograms
    histDataD0 = listOfPrimitivesD0.At(2)

    # PDFs
    pdfD0 = listOfPrimitivesD0.At(3)
    pdfD0S1S2 = listOfPrimitivesD0.At(4)
    pdfD0S1B2 = listOfPrimitivesD0.At(5)
    pdfD0B1S2 = listOfPrimitivesD0.At(6)
    pdfD0B1B2 = listOfPrimitivesD0.At(7)
    pdfD0S1F2 = listOfPrimitivesD0.At(8)
    pdfD0B1F2 = listOfPrimitivesD0.At(9)
    if config["fit"]["add_psi2s"]:
        pdfD0S2S3 = listOfPrimitivesD0.At(10)
        pdfD0B2S3 = listOfPrimitivesD0.At(11)
        pdfD0F2S3 = listOfPrimitivesD0.At(12)

    # pdfD0S1S2.SetLineColor(ROOT.kRed+1)
    # pdfD0S1B2.SetLineColor(ROOT.kAzure+1)
    # pdfD0B1S2.SetLineColor(ROOT.kGreen+1)
    # pdfD0B1B2.SetLineColor(ROOT.kOrange+7)
    
    pdfD0S1S2.SetLineStyle(5)
    pdfD0S1B2.SetLineStyle(5)
    pdfD0B1S2.SetLineStyle(2)
    pdfD0B1B2.SetLineStyle(2)
    pdfD0S1F2.SetLineStyle(5)
    pdfD0B1F2.SetLineStyle(2)
    if config["fit"]["add_psi2s"]:
        pdfD0S2S3.SetLineStyle(3)
        pdfD0B2S3.SetLineStyle(3)
        pdfD0F2S3.SetLineStyle(3)

    canvasOutD0 = ROOT.TCanvas("canvasOutD0", "canvasOutD0", 800, 800)
    canvasOutD0.SetTickx(1)
    canvasOutD0.SetTicky(1)
    frameD0.Draw("AXIS")
    histDataD0.Draw("EP SAME")
    pdfD0.Draw("SAME")
    pdfD0S1S2.Draw("SAME")
    pdfD0S1B2.Draw("SAME")
    pdfD0B1S2.Draw("SAME")
    pdfD0B1B2.Draw("SAME")
    pdfD0S1F2.Draw("SAME")
    pdfD0B1F2.Draw("SAME")
    if config["fit"]["add_psi2s"]:
        pdfD0S2S3.Draw("SAME")
        pdfD0B2S3.Draw("SAME")
        pdfD0F2S3.Draw("SAME")

        
    legend1 = ROOT.TLegend(0.65, 0.65, 0.89, 0.94, " ", "brNDC")
    if config["fit"]["add_psi2s"]:
        legend1 = ROOT.TLegend(0.6, 0.55, 0.8, 0.95, " ", "brNDC")
        
    SetLegend(legend1)
    legend1.SetTextSize(0.030)
    legend1.AddEntry(histDataD0, "data", "EP")
    legend1.AddEntry(pdfD0, "total fit", "L")
    legend1.AddEntry(pdfD0S1S2, "sig. J/#psi - sig. D^{0}", "L")
    legend1.AddEntry(pdfD0S1B2, "sig. J/#psi - bkg. #piK", "L")
    legend1.AddEntry(pdfD0S1F2, "sig. J/#psi - refl. #piK", "L")
    if config["fit"]["add_psi2s"]:
        legend1.AddEntry(pdfD0S2S3, "sig. #psi(2S) - sig. D^{0}", "L")
        legend1.AddEntry(pdfD0B2S3, "sig. #psi(2S) - bkg. #piK", "L")
        legend1.AddEntry(pdfD0F2S3, "sig. #psi(2S) - refl. #piK", "L")
        legend1.AddEntry(pdfD0B1S2, "bkg. "+decayLabel+" - sig. D^{0}", "L")
        legend1.AddEntry(pdfD0B1B2, "bkg. "+decayLabel+" - bkg. #piK", "L")
        legend1.AddEntry(pdfD0B1F2, "bkg. "+decayLabel+" - refl. #piK", "L")
    else:
        legend1.AddEntry(pdfD0B1S2, "bkg. "+decayLabel+" - sig. D^{0}", "L")
        legend1.AddEntry(pdfD0B1B2, "bkg. "+decayLabel+" - bkg. #piK", "L")
        legend1.AddEntry(pdfD0B1F2, "bkg. "+decayLabel+" - refl. #piK", "L")

    legend1.Draw()
    
    latexTitle.DrawLatex(0.18, 0.89, "Work in progress")
    latexTitle.DrawLatex(0.18, 0.83, "pp, #sqrt{#it{s}} = 13.6 TeV")
    latexDecay.DrawLatex(0.18, 0.77, "D^{0} #rightarrow #pi^{+}K^{#minus}")
    if config["fit"]["JpsiChannel"] == "Jpsi2mumu":
        latexDecay.DrawLatex(0.18, 0.73, "J/#psi #rightarrow #mu^{+}#mu^{#minus}")
    else:
        latexDecay.DrawLatex(0.18, 0.73, "J/#psi #rightarrow e^{+}e^{#minus}")
    
    latexRap.DrawLatex(0.18, 0.68, "|#it{#eta}_{#piK}| < 0.8")
    if config["fit"]["JpsiChannel"] == "Jpsi2mumu":
        latexRap.DrawLatex(0.18, 0.63, f"{config['fit']['min_jpsi_rap']} < #it{{#eta}}_{{#mu#mu}} < {config['fit']['max_jpsi_rap']}") #titleSuffix
    else:
        latexRap.DrawLatex(0.18, 0.63, f"|#it{{#eta}}_{{ee}}| < {config['fit']['max_jpsi_rap']}")

    canvasOutD0.Update()
    canvasOutD0.SaveAs(f'{config["output"]["figures"]}/fit_{config["fit"]["JpsiChannel"]}_d0_projection.png')
    
    ## plot the 2D results
    hist3D = listOfPrimitives3D.At(0)
    if config["fit"]["JpsiChannel"] == "Jpsi2mumu":
        hist3D.SetTitle(f';#it{{m}}_{{#piK}} (GeV/#it{{c}}^{{2}});#it{{m}}_{{#mu#mu}} (GeV/#it{{c}}^{{2}})')
        
    else:
        hist3D.SetTitle(f';#it{{m}}_{{#piK}} (GeV/#it{{c}}^{{2}});#it{{m}}_{{ee}} (GeV/#it{{c}}^{{2}})')

    hist3D.GetZaxis().SetTitle(f'Counts per ({round(1000*((config["fit"]["max_fit_range_d0"]-config["fit"]["min_fit_range_d0"])/config["plot_results"]["dataBins"]))} MeV/#it{{c}}^{{2}}#times {round(1000*((config["fit"]["max_fit_range_jpsi"]-config["fit"]["min_fit_range_jpsi"])/config["plot_results"]["dataBins"]))} MeV/#it{{c}}^{{2}})')
    model = listOfPrimitives3D.At(1)
    
    canvasOut3D = ROOT.TCanvas("canvasOut3D", "canvasOut3D", 1200, 1000)
    # canvasOut3D.SetMargin(0.15, 0.15, 0.15, 0.15)
    canvasOut3D.SetTopMargin(0.15)
    ROOT.gStyle.SetOptStat(0)
    hist3D.Draw("LEGO2")
    model.Draw("SURF SAME")
    
    canvasOut3D.SetPhi(230)
    canvasOut3D.SetTheta(20)
    
    latexTitle.DrawLatex(0.1, 0.94, "Work in progress")
    latexTitle.DrawLatex(0.1, 0.88, "pp, #sqrt{#it{s}} = 13.6 TeV")

    latexDecay.DrawLatex(0.1, 0.82, "D^{0} #rightarrow #pi^{+}K^{#minus}")
    if config["fit"]["JpsiChannel"] == "Jpsi2mumu":
        latexDecay.DrawLatex(0.1, 0.78, "J/#psi #rightarrow #mu^{+}#mu^{#minus}")
    else:
        latexDecay.DrawLatex(0.1, 0.78, "J/#psi #rightarrow e^{+}e^{#minus}")

    latexRap.SetTextSize(0.04)
    latexRap.DrawLatex(0.7, 0.92, "|#it{#eta}_{#piK}| < 0.8")
    if config["fit"]["JpsiChannel"] == "Jpsi2mumu":
        latexRap.DrawLatex(0.7, 0.86, f"{config['fit']['min_jpsi_rap']} < #it{{#eta}}_{{#mu#mu}} < {config['fit']['max_jpsi_rap']}") #titleSuffix
    else:
        latexRap.DrawLatex(0.7, 0.86, f"|#it{{#eta}}_{{ee}}| < {config['fit']['max_jpsi_rap']}")

    canvasOut3D.Update()
    canvasOut3D.SaveAs(f'{config["output"]["figures"]}/fit_2d_{config["fit"]["JpsiChannel"]}.png')
    
    input()

if __name__ == '__main__':
    main()

