import sys
import argparse
import yaml
import ROOT
import numpy as np
from array import array
import os

import math
from statistics import mean

sys.path.append('../utils')
sys.path.append('./')

from plot_library import LoadStyle, SetGraStat, SetGraSyst, SetLegend

from model_and_label import getLabel, getGlobalLabel, constructWorkspace

def multi_trial(config, tailParamSet):
    """
    function to fit the 2D distribution of D0 - J/psi masses doing trial variation
    """    
    if config["fit"]["JpsiChannel"] == "Jpsi2ee":
        titleSuffix = "e^{+}e^{-}"
    elif config["fit"]["JpsiChannel"] == "Jpsi2mumu":
        titleSuffix = "#mu^{+}#mu^{-}"
    else:
        print("Error: JpsiChannel not defined in the configuration file.")
        sys.exit(1)

    LoadStyle()
    ROOT.gStyle.SetOptStat(False)
    
    # Variables
    minFitRangeD0 = config["fit"]["min_fit_range_d0"]
    maxFitRangeD0 = config["fit"]["max_fit_range_d0"]
        
    if config["fit"]["unbinned"]:
        fIn = ROOT.TFile(config["inputs"]["data"], "READ")
        sample = fIn.Get(config["inputs"]["tree"])
    else:
        fIn = ROOT.TFile(config["inputs"]["data"], "READ")
        sample = fIn.Get(config["inputs"]["hist"])


    nEvent = sample.GetEntries()
    if config['fit']['weighted'] and config["fit"]["unbinned"]:
        h = ROOT.TH1F("h", "", 1, 0, 1e5)
        sample.Draw("1>>h", "weight")
        nEvent = h.Integral()

    minVarFitRangeD0 = config["fit"]["var_min_fit_range_d0"]
    maxVarFitRangeD0 = config["fit"]["var_max_fit_range_d0"]
    minVarFitRangeJpsi = config["fit"]["var_min_fit_range_jpsi"]
    maxVarFitRangeJpsi = config["fit"]["var_max_fit_range_jpsi"]

    varD0ReflRatio = config["fit"]["var_d0_refl_ratio"]

    for iVarD0RelRatio, d0ReflRatio in enumerate(varD0ReflRatio):
        print(f'iVarD0RelRatio: {iVarD0RelRatio} ; d0ReflRatio: {d0ReflRatio}')

        nTrials = len(minVarFitRangeD0) * len(minVarFitRangeJpsi)
        histMultiTrial = ROOT.TH1D("histMultiTrial", ";Trial;N_{J/#psi-D^{0}}", nTrials, 0, nTrials)
        histMultiTrial.SetMarkerStyle(20)
        histMultiTrial.SetMarkerColor(ROOT.kBlack)
        histMultiTrial.SetLineColor(ROOT.kBlack)

        histCovMatrStatus = ROOT.TH1D("histCovMatrStatus", ";Trial;N_{J/#psi-D^{0}}", nTrials, 0, nTrials)
        histCovMatrStatus.SetMarkerStyle(20)
        histCovMatrStatus.SetMarkerColor(ROOT.kBlack)
        histCovMatrStatus.SetLineColor(ROOT.kBlack)


        # Build fOut name
        directory = config["output"]["directory"]
        
        vecPdfs = ["cb_par_jpsi", "cb_par_psi2s", "cheby_par_jpsi", "cb_par_d0", "cheby_par_d0"]
        strFreePars = config["fit"]["pdf_jpsi_name"] + "_" + config["fit"]["bkg_jpsi_name"] + "_" + config["fit"]["pdf_d0_name"] + "_" + config["fit"]["bkg_d0_name"] + "_d0_refl_ratio_" + str(d0ReflRatio)

        for pdf in vecPdfs:
            vecFreePars = []
            for iPar,par in enumerate(config["fit"][f"{pdf}_is_const"]):
                if not par:
                    vecFreePars.append(iPar)
            if len(vecFreePars) > 0:
                strFreePars = strFreePars + "_" + pdf
                for par in vecFreePars:
                    strFreePars = strFreePars + "_" + str(par)
                strFreePars = strFreePars + "_free"

        fOutName = f'{directory}/multi_trial_{strFreePars}_{tailParamSet}.root'

        fOut = ROOT.TFile(fOutName, "RECREATE")
        iTrial = 0
        parValArray = []
        trialIndexArray  = array( 'f', [] )

        for iVarRangeD0, (minVarRangeD0, maxVarRangeD0) in enumerate(zip(minVarFitRangeD0, maxVarFitRangeD0)):
            print(f'iVarRangeD0: {iVarRangeD0} ; minVarRangeD0: {minVarRangeD0} ; maxVarRangeD0: {maxVarRangeD0}')
            for iVarRangeJpsi, (minVarRangeJpsi, maxVarRangeJpsi) in enumerate(zip(minVarFitRangeJpsi, maxVarFitRangeJpsi)):
                print(f'iVarRangeJpsi: {iVarRangeJpsi} ; minVarRangeJpsi: {minVarRangeJpsi} ; maxVarRangeJpsi: {maxVarRangeJpsi}')
                varArray = [iVarRangeD0, iVarRangeJpsi, iVarD0RelRatio]
                print(iTrial, ") Variation array = ", varArray)

                workSpace = constructWorkspace(config, includeJpsi=True, includeD0=True, nEvent=nEvent, variations=varArray)
                model =  workSpace.pdf("model")

                if config["fit"]["unbinned"]:
                    if config["fit"]["weighted"]:
                        sampleToFit = ROOT.RooDataSet("dataTree", "dataset with weights", sample, ROOT.RooArgSet(workSpace.var("fMassDmes"), workSpace.var("fMass"), workSpace.var("fPtDmes"),  workSpace.var("fPtJpsi"), workSpace.var("fDeltaY"), workSpace.var("fRapJpsi"), workSpace.var("weight")), "", "weight")
                    else:
                        sampleToFit = ROOT.RooDataSet("dataTree", "dataTree", [workSpace.var("fMassDmes"), workSpace.var("fMass"), workSpace.var("fPtDmes"),  workSpace.var("fPtJpsi"), workSpace.var("fDeltaY"), workSpace.var("fRapJpsi")], Import=sample)
                else:
                    sampleToFit = ROOT.RooDataHist("dataHist", "dataHist", [workSpace.var("fMassDmes"), workSpace.var("fMass"), workSpace.var("fPtDmes"),  workSpace.var("fPtJpsi"), workSpace.var("fDeltaY"), workSpace.var("fRapJpsi")], Import=sample)

                if config["fit"]["unbinned"]:
                    ## jpsi rapidity cut
                    sampleToFit = sampleToFit.reduce(f"fRapJpsi>{config['fit']['min_jpsi_rap']} && fRapJpsi<{config['fit']['max_jpsi_rap']}")
                    ## D0 pt cut
                    sampleToFit = sampleToFit.reduce(f"fPtDmes>{config['fit']['min_d0_pt']}")
                    ## dRap cut
                    sampleToFit = sampleToFit.reduce(f"fabs(fDeltaY)>{config['fit']['min_dRap']} && fabs(fDeltaY)<{config['fit']['max_dRap']}")
                    
                # Fit
                fitResult = model.fitTo(sampleToFit, 
                                        ROOT.RooFit.PrintLevel(3), 
                                        ROOT.RooFit.Optimize(1), 
                                        ROOT.RooFit.Hesse(1), 
                                        ROOT.RooFit.Strategy(2), 
                                        ROOT.RooFit.Save(1), 
                                        ROOT.RooFit.Minos(not config['fit']['weighted']), 
                                        ROOT.RooFit.SumW2Error(config['fit']['weighted']), 
                                        ROOT.RooFit.PrintLevel(-1),
                                        ROOT.RooFit.PrintEvalErrors(-1))

                nJPsiD0 = fitResult.floatParsFinal().find("nJPsiD0")
                nBkgJPsi = fitResult.floatParsFinal().find("nBkgJPsi")
                nBkgD0 = fitResult.floatParsFinal().find("nBkgD0")
                nBkgBkg = fitResult.floatParsFinal().find("nBkgBkg")
                print(f'n J/psi - D0 = {nJPsiD0.getVal():.0f} +/- {nJPsiD0.getError():.0f}')
                print(f'Cov. matrix status: {fitResult.covQual()}')
                fitResult.Print()


                # Drawing
                modelHist = model.createHistogram(f'model m_D0, m_J/psi, Trial {iTrial}', workSpace.var("fMassDmes"), Binning=config["plot_results"]["dataBins"], YVar=dict(var=workSpace.var("fMass"), Binning=config["plot_results"]["dataBins"]))
                modelHist.SetLineColor(ROOT.kRed)
                modelHist.SetLineWidth(1)
                if config["fit"]["JpsiChannel"] == "Jpsi2mumu":
                    modelHist.SetTitle(f';#it{{m}}_{{#piK}} (GeV/#it{{c}}^{{2}});#it{{m}}_{{#mu#mu}} (GeV/#it{{c}}^{{2}})')
                else:
                    modelHist.SetTitle(f';#it{{m}}_{{#piK}} (GeV/#it{{c}}^{{2}});#it{{m}}_{{ee}} (GeV/#it{{c}}^{{2}})')

                mJpsiframe = workSpace.var("fMass").frame(Title=" ")
                sampleToFit.plotOn(mJpsiframe, ROOT.RooFit.Name("data"), ROOT.RooFit.Binning(config["plot_results"]["dataBins"]))
                model.plotOn(mJpsiframe, Name={"model"})
                model.plotOn(mJpsiframe, Name={"sigD0_sigJpsiPdf"}, Components={"sigD0_sigJpsiPdf"}, LineStyle=5, LineColor=ROOT.kRed+1)
                model.plotOn(mJpsiframe, Name={"bkgD0_sigJpsiPdf"}, Components={"bkgD0_sigJpsiPdf"}, LineStyle=5, LineColor=ROOT.kAzure+1)
                model.plotOn(mJpsiframe, Name={"bkgJpsi_sigD0Pdf"}, Components={"bkgJpsi_sigD0Pdf"}, LineStyle=2, LineColor=ROOT.kGreen+1)
                model.plotOn(mJpsiframe, Name={"bkgD0_bkgJpsiPdf"}, Components={"bkgD0_bkgJpsiPdf"}, LineStyle=2, LineColor=ROOT.kMagenta+1)
                model.plotOn(mJpsiframe, Name={"reflD0_sigJpsiPdf"}, Components={"reflD0_sigJpsiPdf"}, LineStyle=5, LineColor=ROOT.kYellow+2)
                model.plotOn(mJpsiframe, Name={"reflD0_bkgJpsiPdf"}, Components={"reflD0_bkgJpsiPdf"}, LineStyle=2, LineColor=ROOT.kGray+1)
                if config["fit"]["add_psi2s"]:
                    model.plotOn(mJpsiframe, Name={"sigD0_sigPsi2SPdf"}, Components={"sigD0_sigPsi2SPdf"}, LineStyle=3, LineColor=ROOT.kRed+1)
                    model.plotOn(mJpsiframe, Name={"bkgD0_sigPsi2SPdf"}, Components={"bkgD0_sigPsi2SPdf"}, LineStyle=3, LineColor=ROOT.kAzure+1)
                    model.plotOn(mJpsiframe, Name={"reflD0_sigPsi2SPdf"}, Components={"reflD0_sigPsi2SPdf"}, LineStyle=3, LineColor=ROOT.kYellow+2)

                mD0frame = workSpace.var("fMassDmes").frame(Title=" ")
                sampleToFit.plotOn(mD0frame, ROOT.RooFit.Name("data"), ROOT.RooFit.Binning(config["plot_results"]["dataBins"]))
                model.plotOn(mD0frame, Name={"model"})
                model.plotOn(mD0frame, Name={"sigD0_sigJpsiPdf"}, Components={"sigD0_sigJpsiPdf"}, LineStyle=5, LineColor=ROOT.kRed+1)
                model.plotOn(mD0frame, Name={"bkgD0_sigJpsiPdf"}, Components={"bkgD0_sigJpsiPdf"}, LineStyle=5, LineColor=ROOT.kAzure+1)
                model.plotOn(mD0frame, Name={"bkgJpsi_sigD0Pdf"}, Components={"bkgJpsi_sigD0Pdf"}, LineStyle=2, LineColor=ROOT.kGreen+1)
                model.plotOn(mD0frame, Name={"bkgD0_bkgJpsiPdf"}, Components={"bkgD0_bkgJpsiPdf"}, LineStyle=2, LineColor=ROOT.kMagenta+1)
                model.plotOn(mD0frame, Name={"reflD0_sigJpsiPdf"}, Components={"reflD0_sigJpsiPdf"}, LineStyle=5, LineColor=ROOT.kYellow+2)
                model.plotOn(mD0frame, Name={"reflD0_bkgJpsiPdf"}, Components={"reflD0_bkgJpsiPdf"}, LineStyle=2, LineColor=ROOT.kGray+1)
                if config["fit"]["add_psi2s"]:
                    model.plotOn(mD0frame, Name={"sigD0_sigPsi2SPdf"}, Components={"sigD0_sigPsi2SPdf"}, LineStyle=3, LineColor=ROOT.kRed+1)
                    model.plotOn(mD0frame, Name={"bkgD0_sigPsi2SPdf"}, Components={"bkgD0_sigPsi2SPdf"}, LineStyle=3, LineColor=ROOT.kAzure+1)
                    model.plotOn(mD0frame, Name={"reflD0_sigPsi2SPdf"}, Components={"reflD0_sigPsi2SPdf"}, LineStyle=3, LineColor=ROOT.kYellow+2)
                
                legend_comp = ROOT.TLegend(0.59, 0.57, 0.89, 0.87)
                legend_comp.SetBorderSize(0)
                legend_comp.SetFillStyle(0)
                legend_comp.AddEntry(mD0frame.findObject("data"), "data", "PE")
                legend_comp.AddEntry(mD0frame.findObject("model"), "total fit", "L")

                if config["fit"]["add_psi2s"]:
                    hist_Jpsi = ROOT.TH1F(f'hist_Jpsi_trial_{iTrial}',"", 100, minFitRangeD0, maxFitRangeD0); hist_Jpsi.SetLineWidth(3); hist_Jpsi.SetLineStyle(5); hist_Jpsi.SetLineColor(ROOT.kBlack)
                    hist_psi2s = ROOT.TH1F(f'hist_psi2s_{iTrial}',"", 100, minFitRangeD0, maxFitRangeD0); hist_psi2s.SetLineWidth(3); hist_psi2s.SetLineStyle(3); hist_psi2s.SetLineColor(ROOT.kBlack)
                    hist_sigD0_sigJpsi = ROOT.TH1F(f'hist_sigD0_sigJpsi_{iTrial}',"", 100, minFitRangeD0, maxFitRangeD0); hist_sigD0_sigJpsi.SetLineWidth(3); hist_sigD0_sigJpsi.SetLineStyle(5); hist_sigD0_sigJpsi.SetLineColor(ROOT.kRed+1)
                    hist_bkgD0_sigJpsi = ROOT.TH1F(f'hist_bkgD0_sigJpsi_{iTrial}',"", 100, minFitRangeD0, maxFitRangeD0); hist_bkgD0_sigJpsi.SetLineWidth(3); hist_bkgD0_sigJpsi.SetLineStyle(5); hist_bkgD0_sigJpsi.SetLineColor(ROOT.kAzure+1)
                    hist_bkgJpsi_sigD0 = ROOT.TH1F(f'hist_bkgJpsi_sigD0_{iTrial}',"", 100, minFitRangeD0, maxFitRangeD0); hist_bkgJpsi_sigD0.SetLineWidth(3); hist_bkgJpsi_sigD0.SetLineStyle(2); hist_bkgJpsi_sigD0.SetLineColor(ROOT.kGreen+1)
                    hist_bkgD0_bkgJpsi = ROOT.TH1F(f'hist_bkgD0_bkgJpsi_{iTrial}',"", 100, minFitRangeD0, maxFitRangeD0); hist_bkgD0_bkgJpsi.SetLineWidth(3); hist_bkgD0_bkgJpsi.SetLineStyle(2); hist_bkgD0_bkgJpsi.SetLineColor(ROOT.kMagenta+1)
                    hist_reflD0_sigJpsi = ROOT.TH1F(f'hist_reflD0_sigJpsi_{iTrial}',"", 100, minFitRangeD0, maxFitRangeD0); hist_reflD0_sigJpsi.SetLineWidth(3); hist_reflD0_sigJpsi.SetLineStyle(5); hist_reflD0_sigJpsi.SetLineColor(ROOT.kYellow+2)
                    hist_reflD0_bkgJpsi = ROOT.TH1F(f'hist_reflD0_bkgJpsi_{iTrial}',"", 100, minFitRangeD0, maxFitRangeD0); hist_reflD0_bkgJpsi.SetLineWidth(3); hist_reflD0_bkgJpsi.SetLineStyle(2); hist_reflD0_bkgJpsi.SetLineColor(ROOT.kGray+1)
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

                latexTitle = ROOT.TLatex()
                latexTitle.SetTextSize(0.04)
                latexTitle.SetNDC()
                latexTitle.SetTextFont(42)

                # Fit Jpsi
                canvasFit = ROOT.TCanvas(f'Trial_{iTrial}', f'Trial_{iTrial}', 1600, 800)
                canvasFit.SetTickx(1)
                canvasFit.SetTicky(1)
                canvasFit.Divide(2, 1)

                #mJpsiframe.GetYaxis().SetRangeUser(config["plot_results"]["jpsiFrame"]["y_range"][0], config["plot_results"]["jpsiFrame"]["y_range"][1])

                canvasFit.cd(1)
                mJpsiframe.GetYaxis().SetRangeUser(0, 0.1*sampleToFit.sumEntries() if config['fit']['weighted'] else 0.1*sampleToFit.numEntries() )
                mJpsiframe.GetYaxis().SetTitleOffset(1.4)
                mJpsiframe.GetYaxis().SetTitle(f'Counts per {round(1000*((config["fit"]["max_fit_range_jpsi"]-config["fit"]["min_fit_range_jpsi"])/config["plot_results"]["dataBins"]))} MeV/#it{{c}}^{{2}}')
                mJpsiframe.Draw()
                legend_comp.Draw()
                if config["fit"]["add_psi2s"]:
                    legend_comp2.Draw()

                # Fit D0
                canvasFit.cd(2)
                mD0frame.GetYaxis().SetTitleOffset(1.4)
                mD0frame.GetYaxis().SetLabelSize(0.03)
                mD0frame.GetYaxis().SetRangeUser(0, 0.1*sampleToFit.sumEntries() if config['fit']['weighted'] else 0.1*sampleToFit.numEntries())

                mD0frame.GetYaxis().SetTitle(f'Counts per {round(1000*((config["fit"]["max_fit_range_d0"]-config["fit"]["min_fit_range_d0"])/config["plot_results"]["dataBins"]))} MeV/#it{{c}}^{{2}}')
                mD0frame.Draw()

                latexTitle.DrawLatex(0.3, 0.85, f'#it{{N}}_{{sig. J/#psi - sig. D^{{0}}}} = {nJPsiD0.getVal():.0f} #pm {nJPsiD0.getError():.0f}')
                latexTitle.DrawLatex(0.3, 0.80, f'#it{{N}}_{{sig. J/#psi - bkg. D^{{0}}}} = {nBkgJPsi.getVal():.0f} #pm {nBkgJPsi.getError():.0f}')
                latexTitle.DrawLatex(0.3, 0.75, f'#it{{N}}_{{bkg. J/#psi - sig. D^{{0}}}} = {nBkgD0.getVal():.0f} #pm {nBkgD0.getError():.0f}')
                latexTitle.DrawLatex(0.3, 0.70, f'#it{{N}}_{{bkg. J/#psi - bkg. D^{{0}}}} = {nBkgBkg.getVal():.0f} #pm {nBkgBkg.getError():.0f}')

                canvasFit.Update()
                
                fOut.cd()
                canvasFit.Write()
                
                parValArray.append(nJPsiD0.getVal())
                histMultiTrial.SetBinContent(iTrial+1, nJPsiD0.getVal())
                histMultiTrial.SetBinError(iTrial+1, nJPsiD0.getError())

                histCovMatrStatus.SetBinContent(iTrial+1, fitResult.covQual())
                histCovMatrStatus.SetBinError(iTrial+1, 0)


                iTrial = iTrial + 1
                trialIndexArray.append(iTrial)

                del canvasFit

        #centralVal = mean(parValArray)

        #trialIndexWidthArray = array( 'f', [] )
        #parValSystArray = array( 'f', [] )
        #parErrSystArray = array( 'f', [] )
        #for i in range(0, len(parValArray)):
            #trialIndexWidthArray.append(0.5)
            #parValSystArray.append(mean(parValArray))
            #parErrSystArray.append(ComputeRMS(parValArray))

        #graParSyst = ROOT.TGraphErrors(len(parValArray), trialIndexArray, parValSystArray, trialIndexWidthArray, parErrSystArray)
        #graParSyst.SetFillColorAlpha(ROOT.kGray+1, 0.3)

        #linePar = ROOT.TLine(0, centralVal, len(trialIndexArray), centralVal)
        #linePar.SetLineColor(ROOT.kRed)
        #linePar.SetLineWidth(2)

        #canvasMultiTrial = ROOT.TCanvas("canvasMultiTrial", "", 800, 600)
        #histMultiTrial.GetYaxis().SetRangeUser(0, 1e5)
        #histMultiTrial.Draw("EP")
        #linePar.Draw("SAME")
        #graParSyst.Draw("E2 SAME")
        #canvasMultiTrial.Update()

        #input()

        print(f'[INFO] Writing output on {fOutName} ...')
        histMultiTrial.Write()
        histCovMatrStatus.Write()
        fOut.Close()

    input()
    exit()

def combine_systematics(fileList, dOutName):
    LoadStyle()
    ROOT.gStyle.SetOptStat(False)

    vecYield = []
    vecStatErrYield = []
    nTrialGroups = 0
    histTrialPerGroup = []

    with open(fileList, "r") as fInList:
        nBinsX = 0
        for iFinName,fInName in enumerate(fInList):
            nTrialGroups = nTrialGroups + 1
            fIn = ROOT.TFile(fInName.strip(), "READ")
            histMultiTrial = fIn.Get("histMultiTrial")
            histCovMatrStatus = fIn.Get("histCovMatrStatus")

            offset = nBinsX
            nBinsX = nBinsX + histMultiTrial.GetNbinsX()

            histTrialPerGroup.append(ROOT.TH1D(f'histMultiTrialGroup{iFinName}', ";Trial;N_{J/#psi-D^{0}}", nBinsX, 0, nBinsX))
            histTrialPerGroup[iFinName].SetDirectory(0)
            for iBin in range(0, nBinsX):
                rawYield = float(histMultiTrial.GetBinContent(iBin+1))
                errRawYield = float(histMultiTrial.GetBinError(iBin+1))
                if rawYield == 0: continue 
                relErr = errRawYield / rawYield

                if histCovMatrStatus.GetBinContent(iBin+1) == 3 and relErr < 0.5:
                    vecYield.append(histMultiTrial.GetBinContent(iBin+1))
                    vecStatErrYield.append(histMultiTrial.GetBinError(iBin+1))
                    histTrialPerGroup[iFinName].SetBinContent(iBin+1+offset, histMultiTrial.GetBinContent(iBin+1))
                    histTrialPerGroup[iFinName].SetBinError(iBin+1+offset, histMultiTrial.GetBinError(iBin+1))

    for i in range(nTrialGroups):
        col_index = int( (i / max(1, nTrialGroups-1)) * 255 )
        histTrialPerGroup[i].SetMarkerColor(ROOT.TColor.GetColorPalette(col_index))
        histTrialPerGroup[i].SetLineColor(ROOT.TColor.GetColorPalette(col_index))

    histYield = ROOT.TH1D("histYield", ";N_{J/#psi-D^{0}}", 100, 0, 2 * mean(vecYield))
    histYield.SetLineColor(ROOT.kAzure+4)
    histYield.SetFillStyle(0)
    histYield.SetFillColorAlpha(ROOT.kAzure+4, 0.5)

    histStatErrYield = ROOT.TH1D("histStatErrYield", ";(stat. error)_{J/#psi-D^{0}}", 100, 0, 2 * mean(vecStatErrYield))
    histStatErrYield.SetLineColor(ROOT.kAzure+4)
    histStatErrYield.SetFillStyle(0)
    histStatErrYield.SetFillColorAlpha(ROOT.kAzure+4, 0.5)

    for (val, err) in zip(vecYield, vecStatErrYield):
        histYield.Fill(val)
        histStatErrYield.Fill(err)

    latexTitle = ROOT.TLatex()
    latexTitle.SetTextSize(0.04)
    latexTitle.SetNDC()
    latexTitle.SetTextFont(42)

    meanYield = histYield.GetMean()
    meanStatErrYield = histStatErrYield.GetMean()
    RmsYield = histYield.GetRMS()

    lineYieldMean = ROOT.TLine(meanYield, 0, meanYield, histYield.GetMaximum())
    lineYieldMean.SetLineStyle(ROOT.kSolid)
    lineYieldMean.SetLineColor(ROOT.kRed+1)
    lineYieldMean.SetLineWidth(2)

    lineYieldRmsUp = ROOT.TLine(meanYield+RmsYield, 0, meanYield+RmsYield, histYield.GetMaximum())
    lineYieldRmsUp.SetLineStyle(ROOT.kDashed)
    lineYieldRmsUp.SetLineColor(ROOT.kRed+1)
    lineYieldRmsUp.SetLineWidth(2)

    lineYieldRmsLow = ROOT.TLine(meanYield-RmsYield, 0, meanYield-RmsYield, histYield.GetMaximum())
    lineYieldRmsLow.SetLineStyle(ROOT.kDashed)
    lineYieldRmsLow.SetLineColor(ROOT.kRed+1)
    lineYieldRmsLow.SetLineWidth(2)
    
    canvasYield = ROOT.TCanvas("canvasYield", "", 800, 600)
    histYield.Draw("H")
    lineYieldMean.Draw("SAME")
    lineYieldRmsUp.Draw("SAME")
    lineYieldRmsLow.Draw("SAME")
    latexTitle.DrawLatex(0.70, 0.85, f'N. trials = {len(vecYield)}')
    latexTitle.DrawLatex(0.70, 0.80, f'Mean = {meanYield:.0f}')
    latexTitle.DrawLatex(0.70, 0.75, f'RMS = {RmsYield:.0f} ({(RmsYield/meanYield)*100:.0f}%)')
    canvasYield.Update()

    canvasStatErrYield = ROOT.TCanvas("canvasStatErrYield", "", 800, 600)
    histStatErrYield.Draw("H")
    latexTitle.DrawLatex(0.70, 0.85, f'N. trials = {len(vecYield)}')
    latexTitle.DrawLatex(0.70, 0.80, f'Mean = {meanStatErrYield:.0f} ({(meanStatErrYield/meanYield)*100:.0f}%)')
    canvasStatErrYield.Update()

    canvasYieldPerGroup = ROOT.TCanvas("canvasYieldPerGroup", "", 800, 600)
    histTrialPerGroup[nTrialGroups-1].GetYaxis().SetRangeUser(0, 2 * meanYield)
    histTrialPerGroup[nTrialGroups-1].SetLineColor(ROOT.kRed)
    histTrialPerGroup[nTrialGroups-1].SetMarkerColor(ROOT.kRed)
    histTrialPerGroup[nTrialGroups-1].Draw("EP")
    for iTrial in range(nTrialGroups-2, -1, -1):
        histTrialPerGroup[iTrial].Draw("EP SAME")
    canvasYieldPerGroup.Update()

    input()

    canvasYield.SaveAs(f'{dOutName}/yield.pdf')
    canvasStatErrYield.SaveAs(f'{dOutName}/statErrYield.pdf')
    canvasYieldPerGroup.SaveAs(f'{dOutName}/yieldPerGroup.pdf')

    exit()

def ComputeRMS(parValArray):
    '''
    Method to evaluate the RMS of a sample ()
    '''
    mean = 0
    for parVal in parValArray:
        mean += parVal
    mean = mean / len(parValArray)
    stdDev = 0
    for parVal in parValArray:
        stdDev += (parVal - mean) * (parVal - mean)
    stdDev = math.sqrt(stdDev / len(parValArray))
    return stdDev



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', nargs='?', default='config.yml', help='config file name')
    parser.add_argument("--run", help="Run multi-trial fit", action="store_true")
    parser.add_argument("--combine", help="Combine multi-trial files [requires filesToMerge]", action="store_true")
    parser.add_argument("-t", "--text", dest="fileList", type=str, help="List of files to be combined", required=False)
    parser.add_argument("-o", "--output", dest="dOutName", type=str, help="Name of the output directory", required=False)
    
    args = parser.parse_args()
    
    if args.run:
        if args.cfgFileName is not None:
            with open(args.cfgFileName, 'r') as yml_cfg:
                inputCfg = yaml.load(yml_cfg, yaml.FullLoader)
            
            tailParamSet = ""
            if "data_tails" in args.cfgFileName:
                tailParamSet = "data_tails"
            if "mc_tails" in args.cfgFileName:
                tailParamSet = "mc_tails"

            multi_trial(inputCfg, tailParamSet)
        else: 
            print(f'[ERROR] No configuration passed to the function')

    if args.combine and args.fileList:
        combine_systematics(args.fileList, args.dOutName)

