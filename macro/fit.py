import argparse
import ROOT
from ROOT import *

def main():
    print('start')
    parser = argparse.ArgumentParser(description='Arguments to pass')
    #parser.add_argument('cfgFileName', metavar='text', default='config.yml', help='config file name')
    parser.add_argument("--do_fit", help="Do the fit to data / toy MC", action="store_true")
    parser.add_argument("--do_prefit", help="Do the pre-fit to data", action="store_true")
    args = parser.parse_args()

    if args.do_fit:
        fit()
    
    if args.do_prefit:
        prefit()


def prefit():
    # Variables
    mD0   = ROOT.RooRealVar("D-meson mass", "#it{m}_{#pi#it{K}} (GeV/#it{c}^{2})", 1.76, 2.00)
    mJpsi = ROOT.RooRealVar("Dimuon mass", "#it{m}_{#mu#mu} (GeV/#it{c}^{2})", 2.50, 4.00)

     # Load data for pre-fit
    fIn = ROOT.TFile("histograms.root", "READ")
    dataHistD0 = ROOT.RooDataHist("dataHistD0", "dataHistD0", [mD0], Import=fIn.Get("hMassDmeson"))
    dataHistJpsi = ROOT.RooDataHist("dataHistJpsi", "dataHistJpsi", [mJpsi], Import=fIn.Get("hMassJPsi"))

    # Pre-fit D0
    meanD0  = ROOT.RooRealVar("meanD0", "meanD0", 1.850, 1.750, 1.950)
    sigmaD0 = ROOT.RooRealVar("sigmaD0", "sigmaD0", 0.050, 0.010, 0.080)
    gausPdfD0 = ROOT.RooGaussian("gausPdfD0", "Crystal Ball D0", mD0, meanD0, sigmaD0)

    chebyParsD0 = [ROOT.RooRealVar(f"cheb_coeff_{i}", f"Coeff_{i}", 0.05, -10, 10) for i in range(3)]
    chebyPdfD0 = ROOT.RooChebychev("chebyPdfD0", "Cheby for Bkg1", mD0, ROOT.RooArgList(*chebyParsD0))

    nSigD0  = ROOT.RooRealVar("nSigD0", "D0 signal", 1e6, 0., 1e8)
    nBkgD0  = ROOT.RooRealVar("nBkgD0", "D0 background", 1e6, 0., 1e8)
    modelD0 = ROOT.RooAddPdf("modelD0", "sigD0 + bkgD0", ROOT.RooArgList(gausPdfD0, chebyPdfD0), ROOT.RooArgList(nSigD0, nBkgD0))

    fitResultD0 = modelD0.fitTo(dataHistD0, ROOT.RooFit.PrintLevel(3), ROOT.RooFit.Optimize(1), ROOT.RooFit.Hesse(1), ROOT.RooFit.Minos(1), ROOT.RooFit.Strategy(2), ROOT.RooFit.Save(1))

    mD0frame = mD0.frame(Title="D-meson")
    dataHistD0.plotOn(mD0frame)
    modelD0.plotOn(mD0frame)
    modelD0.plotOn(mD0frame, Name={"Sig"}, Components={gausPdfD0}, LineStyle="--", LineColor=ROOT.kRed+1)
    modelD0.plotOn(mD0frame, Name={"Bkg"}, Components={chebyPdfD0}, LineStyle="--", LineColor=ROOT.kAzure+4)

    # Pre-fit J/psi
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

    mJpsiframe = mJpsi.frame(Title="Dmuon")
    dataHistJpsi.plotOn(mJpsiframe)
    modelJpsi.plotOn(mJpsiframe)
    modelJpsi.plotOn(mJpsiframe, Name={"Sig"}, Components={cbPdfJpsi}, LineStyle="--", LineColor=ROOT.kRed+1)
    modelJpsi.plotOn(mJpsiframe, Name={"Bkg"}, Components={chebyPdfJpsi}, LineStyle="--", LineColor=ROOT.kAzure+4)

    canvasFit = ROOT.TCanvas("canvasFit", "canvasFit", 1200, 600)
    canvasFit.Divide(2, 1)

    canvasFit.cd(1)
    mD0frame.GetYaxis().SetTitleOffset(1.4)
    mD0frame.Draw()

    canvasFit.cd(2)
    mJpsiframe.GetYaxis().SetTitleOffset(1.4)
    mJpsiframe.Draw()

    canvasFit.Update()

    fitResultD0.Print()
    fitResultJpsi.Print()
    input()
    exit()

def fit():
    # Variables
    mD0   = ROOT.RooRealVar("D-meson mass", "#it{m}_{#pi#it{K}} (GeV/#it{c}^{2})", 1.75, 2.00)
    mJpsi = ROOT.RooRealVar("Dimuon mass", "#it{m}_{#mu#mu} (GeV/#it{c}^{2})", 2.50, 4.00)

    # Yields parameters
    genJpsiD0  = 5000
    genBkgJpsi = 1000
    genBkgD0   = 2000
    genBkgBkg  = 10000

    nJPsiD0  = ROOT.RooRealVar("nJPsiD0", "number of JPsi-D0", genJpsiD0, 0., genJpsiD0 * 5.)
    nBkgJPsi = ROOT.RooRealVar("nBkgJPsi", "number of JPsi-Bkg", genBkgJpsi, 0., genBkgJpsi * 5.)
    nBkgD0   = ROOT.RooRealVar("nBkgD0", "number of D0-Bkg", genBkgD0, 0., genBkgD0 * 5.)
    nBkgBkg  = ROOT.RooRealVar("nBkgBkg", "number of Bkg-Bkg", genBkgBkg, 0., genBkgBkg * 5.)

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

    # Generate toy sample
    data  = sigD0_sigJpsiPdf.generate(ROOT.RooArgSet(mD0, mJpsi), genJpsiD0)
    sig1bkg1Sample  = bkgD0_sigJpsiPdf.generate(ROOT.RooArgSet(mD0, mJpsi), genBkgJpsi)
    sig2bkg2Sample  = bkgJpsi_sigD0Pdf.generate(ROOT.RooArgSet(mJpsi, mD0), genBkgD0)
    bkg1bkg2Sample = bkgD0_bkgJpsiPdf.generate(ROOT.RooArgSet(mD0, mJpsi), genBkgBkg)

    data.append(sig1bkg1Sample)
    data.append(sig2bkg2Sample)
    data.append(bkg1bkg2Sample)

    # Create histograms for plot / fitting
    #histJpsiD0 = data.createHistogram("data m_{J/#psi},m_{D0}", mD0, Binning=50, YVar=dict(var=mJpsi, Binning=50))
    fIn = ROOT.TFile("histograms.root", "READ")
    histJpsiD0= fIn.Get("hSparseJPsiDmeson_proj_3_2")
    dataHist = ROOT.RooDataHist("dataHist", "dataHist", [mD0, mJpsi], Import=histJpsiD0)

    # Fit
    fitResult = model.fitTo(dataHist, ROOT.RooFit.PrintLevel(3), ROOT.RooFit.Optimize(1), ROOT.RooFit.Hesse(1), ROOT.RooFit.Minos(1), ROOT.RooFit.Strategy(2), ROOT.RooFit.Save(1))
    modelHist = model.createHistogram("model m_{D0}, m_{J/#psi}", mD0, Binning=50, YVar=dict(var=mJpsi, Binning=50))
    modelHist.SetLineColor(ROOT.kRed)
    modelHist.SetLineWidth(2)

    mJpsiframe = mJpsi.frame(Title="Dimuon")
    dataHist.plotOn(mJpsiframe)
    model.plotOn(mJpsiframe)
    model.plotOn(mJpsiframe, Name={"sigD0_sigJpsiPdf"}, Components={sigD0_sigJpsiPdf}, LineStyle="--", LineColor=ROOT.kRed+1)
    model.plotOn(mJpsiframe, Name={"bkgD0_sigJpsiPdf"}, Components={bkgD0_sigJpsiPdf}, LineStyle="--", LineColor=ROOT.kAzure+4)
    model.plotOn(mJpsiframe, Name={"bkgJpsi_sigD0Pdf"}, Components={bkgJpsi_sigD0Pdf}, LineStyle="--", LineColor=ROOT.kGreen+1)
    model.plotOn(mJpsiframe, Name={"bkgD0_bkgJpsiPdf"}, Components={bkgD0_bkgJpsiPdf}, LineStyle="--", LineColor=ROOT.kMagenta)

    mD0frame = mD0.frame(Title="D-meson")
    dataHist.plotOn(mD0frame)
    model.plotOn(mD0frame)
    model.plotOn(mD0frame, Name={"sigD0_sigJpsiPdf"}, Components={sigD0_sigJpsiPdf}, LineStyle="--", LineColor=ROOT.kRed+1)
    model.plotOn(mD0frame, Name={"bkgD0_sigJpsiPdf"}, Components={bkgD0_sigJpsiPdf}, LineStyle="--", LineColor=ROOT.kAzure+4)
    model.plotOn(mD0frame, Name={"bkgJpsi_sigD0Pdf"}, Components={bkgJpsi_sigD0Pdf}, LineStyle="--", LineColor=ROOT.kGreen+1)
    model.plotOn(mD0frame, Name={"bkgD0_bkgJpsiPdf"}, Components={bkgD0_bkgJpsiPdf}, LineStyle="--", LineColor=ROOT.kMagenta)

    legend = ROOT.TLegend(0.55, 0.65, 0.85, 0.85)
    legend.AddEntry(mJpsiframe.findObject("sigD0_sigJpsiPdf"), "sig. J/#psi - sig. D0", "L")
    legend.AddEntry(mJpsiframe.findObject("bkgD0_sigJpsiPdf"), "sig. J/#psi - bkg. D0", "L")
    legend.AddEntry(mJpsiframe.findObject("bkgJpsi_sigD0Pdf"), "bkg. J/#psi - sig. D0", "L")
    legend.AddEntry(mJpsiframe.findObject("bkgD0_bkgJpsiPdf"), "bkg. J/#psi - bkg. D0", "L")
    legend.Draw()

    canvasFit = ROOT.TCanvas("canvasFit", "canvasFit", 1200, 600)
    canvasFit.Divide(2, 1)

    canvasFit.cd(1)
    mJpsiframe.GetYaxis().SetTitleOffset(1.4)
    mJpsiframe.Draw()
    legend.Draw()

    canvasFit.cd(2)
    mD0frame.GetYaxis().SetTitleOffset(1.4)
    mD0frame.Draw()
    legend.Draw()
    canvasFit.Update()

    canvasFitHist = ROOT.TCanvas("canvasFitHist", "canvasFitHist", 1000, 1000)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gPad.SetLeftMargin(0.15)
    histJpsiD0.GetXaxis().SetTitleOffset(2.0)
    histJpsiD0.GetYaxis().SetTitleOffset(2.0)
    histJpsiD0.GetZaxis().SetTitleOffset(2.0)
    histJpsiD0.Draw("LEGO2")
    modelHist.Draw("SURF SAME")
    canvasFitHist.Update()

    fitResult.Print()

    fOut = ROOT.TFile("myTest.root", "RECREATE")
    canvasFitHist.Write()
    fOut.Close()

    input()
    exit()

if __name__ == '__main__':
    main()
