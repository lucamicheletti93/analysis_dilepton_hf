void LoadStyle();

void sigma_eff_collection() {
    LoadStyle();

    //------------------------------------------------------------------------------------------------------------------------//
    // Midrapidity results
    double yCentrMidRun3[] = {2.5};
    double yErrMidRun3[] = {0};
    double sigmaEffMidRun3[] = {11.}; // Dummy value, to be raplaced with new one
    double statSigmaEffMidRun3[] = {3.}; // Dummy value, to be raplaced with new one

    TGraphErrors *graSigmaEffMidRun3 = new TGraphErrors(1, sigmaEffMidRun3, yCentrMidRun3, statSigmaEffMidRun3, yErrMidRun3);
    graSigmaEffMidRun3 -> SetMarkerStyle(20);
    graSigmaEffMidRun3 -> SetMarkerColor(kRed+1);
    graSigmaEffMidRun3 -> SetLineColor(kRed+1);

    // Forward rapidity results
    double yCentrFwdRun3[] = {1.5};
    double yErrFwdRun3[] = {0};
    double sigmaEffFwdRun3[] = {11.}; // Dummy value, to be raplaced with new one
    double statSigmaEffFwdRun3[] = {3.}; // Dummy value, to be raplaced with new one

    TGraphErrors *graSigmaEffFwdRun3 = new TGraphErrors(1, sigmaEffFwdRun3, yCentrFwdRun3, statSigmaEffFwdRun3, yErrFwdRun3);
    graSigmaEffFwdRun3 -> SetMarkerStyle(20);
    graSigmaEffFwdRun3 -> SetMarkerColor(kRed+1);
    graSigmaEffFwdRun3 -> SetLineColor(kRed+1);
    //------------------------------------------------------------------------------------------------------------------------//

    double xCentrs[] = {14.5, 12.1, 16.4, 4.8, 2.2, 20.7, 15, 42, 34.8, 8.2, 6.1, 17.2, 6.3, 4.7, 14, 26.3, 13.1, 6.2};
    double yCentrs[] = {25.5, 24.5, 22.5, 21.5, 20.5, 18.5, 17.5, 16.5, 15.5, 14.5, 13.5, 11.5, 10.5, 9.5, 7.5, 6.5, 5.5, 4.5};
    double xErrLow[] = {2.9, 5.4, 2.3, 2.5, 1.14, 6.6, 4.2, 5, 8.4, 2, 1.9, 2.25, 1.9, 1.5, 4.8, 6.4, 2.9, 1.78};
    double xErrHigh[] = {2.4, 10.7, 2.3, 2.5, 6, 6.6, 5.8, 5, 2.6, 2, 3.3, 2.25, 1.9, 2.4, 8.2, 22.8, 2.9, 1.78};
    double yErrLow[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double yErrHigh[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    string yLabels[] = {
        "CDF (#gamma-3 jets)", 
        "CDF (4 jets)",
        "",
        "D0 (#gamma-3 jets)",
        "D0 (J/#psi-J/#psi)",
        "D0 (J/#psi-#varUpsilon)*",
        "",
        "CMS (W^{#pm}-2 jets)",
        "ATLAS (W^{#pm}-2 jets)",
        "LHCb (D^{0}-D^{0})",
        "LHCb (J/#psi-D^{0})*",
        "CMS (J/#psi-J/#psi)*",
        "ATLAS (J/#psi-W^{#pm})*",
        "",
        "LHCb (#varUpsilon(1S)-D^{0})",
        "ATLAS (J/#psi-J/#psi)",
        "ATLAS (J/#psi-Z^{0})*",
        "",
        "LHCb (J/#psi-#varUpsilon(2S))",
        "LHCb (J/#psi-#varUpsilon(1S))",
        "LHCb (J/#psi-J/#psi)",
        "ALICE (J/#psi-J/#psi)",
        "",
        "ALICE (J/#psi-D^{0}, #Deltay < 1.5)",
        "ALICE (J/#psi-D^{0}, 2.1 < #Deltay < 4.6)"
    };

    TGraphAsymmErrors *graSigmaEffToFit = new TGraphAsymmErrors(18, yCentrs, xCentrs, yErrLow, yErrHigh, xErrLow, xErrHigh);
    graSigmaEffToFit -> SetMarkerStyle(24);
    graSigmaEffToFit -> SetMarkerColor(kBlack);
    graSigmaEffToFit -> SetLineColor(kBlack);

    TF1 *funcPol0 = new TF1("funcPol0", "pol0", 0, 30);
    graSigmaEffToFit -> Fit(funcPol0, "R0Q");

    TLine *lineAverage = new TLine(funcPol0 -> GetParameter(0), 0, funcPol0 -> GetParameter(0), 30);
    lineAverage -> SetLineColorAlpha(kAzure+4, 0.6);
    lineAverage -> SetLineWidth(2);

    TCanvas *canvasSigmaEff = new TCanvas("canvasSigmaEff", "canvasSigmaEff", 800, 900);
    canvasSigmaEff -> SetTickx(1);
    canvasSigmaEff -> SetTicky(1);
    canvasSigmaEff -> SetLeftMargin(0.35);
    canvasSigmaEff -> SetRightMargin(0.05);
    canvasSigmaEff -> SetBottomMargin(0.17);
    canvasSigmaEff -> SetFrameLineWidth(2);
    canvasSigmaEff -> SetFrameBorderMode(0);
    canvasSigmaEff -> SetFrameLineWidth(2);
    canvasSigmaEff -> SetFrameBorderMode(0);

    TH2D *histSigmaEff = new TH2D("histSigmaEff", "", 100, 0, 60, 30, 0, 30);
    histSigmaEff -> GetXaxis() -> SetTitle("#sigma_{eff} (mb)");
    histSigmaEff -> GetYaxis() -> SetLabelSize(0.035);
    for (int iBin = 0;iBin < 26;iBin++) {
        histSigmaEff -> GetYaxis() -> SetBinLabel(26-iBin, yLabels[iBin].c_str());
    }
    histSigmaEff -> Draw();

    TGraphAsymmErrors *graSigmaEff = new TGraphAsymmErrors(18, xCentrs, yCentrs, xErrLow, xErrHigh, yErrLow, yErrHigh);
    graSigmaEff -> SetMarkerStyle(20);
    graSigmaEff -> SetMarkerColor(kBlack);
    graSigmaEff -> SetLineColor(kBlack);
    graSigmaEff -> Draw("EP");

    TLine *line0 = new TLine(0, 27, 60, 27);
    line0 -> SetLineColor(kGray+1);
    line0 -> SetLineStyle(kDashed);
    line0 -> Draw();

    TLatex *tex1 = new TLatex(36, 26, "#it{p#bar{p}}, #sqrt{s}=1.8 TeV");
    tex1 -> SetTextFont(42);
    tex1 -> SetTextSize(0.03);
    tex1 -> SetLineWidth(2);
    tex1 -> Draw();

    TLine *line1 = new TLine(0, 24, 60, 24);
    line1 -> SetLineColor(kGray+1);
    line1 -> SetLineStyle(kDashed);
    line1 -> Draw();

    TLatex *tex2 = new TLatex(36, 23, "#it{p#bar{p}}, #sqrt{s}=1.96 TeV");
    tex2 -> SetTextFont(42);
    tex2 -> SetTextSize(0.03);
    tex2 -> SetLineWidth(2);
    tex2 -> Draw();

    TLine *line2 = new TLine(0, 20, 60, 20);
    line2 -> SetLineColor(kGray+1);
    line2 -> SetLineStyle(kDashed);
    line2 -> Draw();

    TLatex *tex3 = new TLatex(36, 19, "#it{pp}, #sqrt{s}=7 TeV");
    tex3 -> SetTextFont(42);
    tex3 -> SetTextSize(0.03);
    tex3 -> SetLineWidth(2);
    tex3 -> Draw();

    TLine *line3 = new TLine(0, 12, 60, 12);
    line3 -> SetLineColor(kGray+1);
    line3 -> SetLineStyle(kDashed);
    line3 -> Draw();

    TLatex *tex4 = new TLatex(36, 11, "#it{pp}, #sqrt{s}=8 TeV");
    tex4 -> SetTextFont(42);
    tex4 -> SetTextSize(0.03);
    tex4 -> SetLineWidth(2);
    tex4 -> Draw();

    TLine *line4 = new TLine(0, 9, 60, 9);
    line4 -> SetLineColor(kGray+1);
    line4 -> SetLineStyle(kDashed);
    line4 -> Draw();

    TLatex *tex5 = new TLatex(36, 8, "#it{pp}, #sqrt{s}=13 TeV");
    tex5 -> SetTextFont(42);
    tex5 -> SetTextSize(0.03);
    tex5 -> SetLineWidth(2);
    tex5 -> Draw();

    TLine *line5 = new TLine(0, 4, 60, 4);
    line5 -> SetLineColor(kGray+1);
    line5 -> SetLineStyle(kDashed);
    line5 -> Draw();

    TLatex *tex6 = new TLatex(36, 3, "#it{pp}, #sqrt{s}=13.6 TeV");
    tex6 -> SetTextFont(42);
    tex6 -> SetTextSize(0.03);
    tex6 -> SetLineWidth(2);
    tex6 -> Draw();

    lineAverage -> Draw();

    TLegend *legend = new TLegend(0.72, 0.89, 0.85, 0.92, " ", "brNDC");
    legend -> SetBorderSize(0);
    legend -> SetFillColor(10);
    legend -> SetFillStyle(1);
    legend -> SetLineStyle(0);
    legend -> SetLineColor(0);
    legend -> SetLineWidth(0);
    legend -> SetTextFont(42);
    legend -> SetTextSize(0.03);
    legend -> AddEntry(lineAverage, "Average", "L");
    legend -> Draw();

    graSigmaEffMidRun3 -> Draw("EP SAME");
    graSigmaEffFwdRun3 -> Draw("EP SAME");

    canvasSigmaEff -> SaveAs("sigma_effective.pdf");
}
////////////////////////////////////////////////////////////////////////////////
void LoadStyle(){
    int font = 42;
    gStyle -> SetFrameBorderMode(0);
    gStyle -> SetFrameFillColor(0);
    gStyle -> SetCanvasBorderMode(0);
    gStyle -> SetPadBorderMode(0);
    gStyle -> SetPadColor(10);
    gStyle -> SetCanvasColor(10);
    gStyle -> SetTitleFillColor(10);
    gStyle -> SetTitleBorderSize(1);
    gStyle -> SetStatColor(10);
    gStyle -> SetStatBorderSize(1);
    gStyle -> SetLegendBorderSize(1);
    gStyle -> SetDrawBorder(0);
    gStyle -> SetTextFont(font);
    gStyle -> SetStatFontSize(0.05);
    gStyle -> SetStatX(0.97);
    gStyle -> SetStatY(0.98);
    gStyle -> SetStatH(0.03);
    gStyle -> SetStatW(0.3);
    gStyle -> SetTickLength(0.02,"y");
    gStyle -> SetEndErrorSize(3);
    gStyle -> SetLabelSize(0.05,"xyz");
    gStyle -> SetLabelFont(font,"xyz");
    gStyle -> SetLabelOffset(0.01,"xyz");
    gStyle -> SetTitleFont(font,"xyz");
    gStyle -> SetTitleOffset(0.9,"x");
    gStyle -> SetTitleOffset(1.02,"y");
    gStyle -> SetTitleSize(0.05,"xyz");
    gStyle -> SetMarkerSize(1.3);
    gStyle -> SetOptStat(0);
    gStyle -> SetEndErrorSize(0);
    gStyle -> SetCanvasPreferGL(kTRUE);
    gStyle -> SetHatchesSpacing(0.5);
    gStyle -> SetPadLeftMargin(0.15);
    gStyle -> SetPadBottomMargin(0.15);
    gStyle -> SetPadTopMargin(0.05);
    gStyle -> SetPadRightMargin(0.05);
    gStyle -> SetEndErrorSize(0.0);
    gStyle -> SetTitleSize(0.05,"X");
    gStyle -> SetTitleSize(0.045,"Y");
    gStyle -> SetLabelSize(0.045,"X");
    gStyle -> SetLabelSize(0.045,"Y");
    gStyle -> SetTitleOffset(1.2,"X");
    gStyle -> SetTitleOffset(1.35,"Y");
}
