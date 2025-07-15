#include <array>
#include <string>
#include <numeric>
#include <map>

#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH3.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TF1.h>
#include <TSpline.h>
#include <TNtuple.h>

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"

#include "Pythia8/Pythia.h"

namespace
{
    enum decayer
    {
        kPythia8 = 0,
        kEvtGen
    };

    enum tunes
    {
        kMonash = 0,
        kCRMode0,
        kCRMode2,
        kCRMode3
    };

    enum processes
    {
        kSoftQCD = 0,
        kHardQCD
    };

    const array<int, 4> absPdgOpenCharm = {411, 421, 431, 4122};
}


template<typename T, typename PPythia>
bool isFromBeauty(T &mothers, PPythia& pythia)
{
    for(auto const& mom : mothers) {
      int absPdgMom = std::abs(pythia.event[mom].id());
      if(absPdgMom == 5 || absPdgMom/100 == 5 || absPdgMom/1000 == 5 ||
        (absPdgMom-10000)/100 == 5 || (absPdgMom-20000)/100 == 5 || (absPdgMom-30000)/100 == 5 ||
        (absPdgMom-100000)/100 == 5 || (absPdgMom-200000)/100 == 5 || (absPdgMom-300000)/100 == 5) {
          return true;
      }
    }

  return false;
}

//__________________________________________________________________________________________________
template<typename PPythia>
void runSimulation(PPythia& pythia, TNtuple* tuplePairs, TH1D* histEvents, std::map<int, TH3F*> histPtVsY, int nEvents)
{
    //__________________________________________________________
    // perform the simulation
    for (auto iEvent{0}; iEvent<nEvents; ++iEvent)
    {
        if(!pythia.next()) {
            continue;
        }

        histEvents->Fill(1.5f);

        std::vector<Pythia8::Particle> openCharms;
        std::vector<Pythia8::Particle> jPsis;
        std::vector<bool> jPsiFromB{};
        std::vector<bool> openCharmFromB{};
        std::vector<std::vector<int>> jPsiMothers{};
        std::vector<std::vector<int>> openCharmMothers{};
        int nFound{0};
        for(auto iPart{2}; iPart<pythia.event.size(); ++iPart)
        {
            int pdg = pythia.event[iPart].id();
            int absPdg = std::abs(pdg);
            if (std::find(absPdgOpenCharm.begin(), absPdgOpenCharm.end(), absPdg) != absPdgOpenCharm.end()) {
                auto mothers = pythia.event[iPart].motherList();
                size_t initMomSize{0u};
                while (mothers.size() > initMomSize) {
                    auto mothersOld = mothers;
                    initMomSize = mothersOld.size();
                    for (auto const& mom : mothersOld) {
                        auto mothersNew = pythia.event[mom].motherList();
                        for (auto const& momNew : mothersNew) {
                            if (std::find(mothers.begin(), mothers.end(), momNew) == mothers.end() && std::abs(pythia.event[momNew].id()) != 2212) {
                                mothers.push_back(momNew);
                            }
                        }
                    }
                }
                bool fromB = isFromBeauty(mothers, pythia);
                histPtVsY[absPdg]->Fill(pythia.event[iPart].pT(), pythia.event[iPart].y(), int(fromB));
                openCharms.push_back(pythia.event[iPart]);
                openCharmFromB.push_back(fromB);
                openCharmMothers.push_back(mothers);
            }  else if (absPdg == 443) {
                auto mothers = pythia.event[iPart].motherList();
                size_t initMomSize{0u};
                while (mothers.size() > initMomSize) {
                    auto mothersOld = mothers;
                    initMomSize = mothersOld.size();
                    for (auto const& mom : mothersOld) {
                        auto mothersNew = pythia.event[mom].motherList();
                        for (auto const& momNew : mothersNew) {
                            if (std::find(mothers.begin(), mothers.end(), momNew) == mothers.end() && std::abs(pythia.event[momNew].id()) != 2212) {
                                mothers.push_back(momNew);
                            }
                        }
                    }
                }
                bool fromB = isFromBeauty(mothers, pythia);
                histPtVsY[absPdg]->Fill(pythia.event[iPart].pT(), pythia.event[iPart].y(), int(fromB));  
                jPsis.push_back(pythia.event[iPart]);
                jPsiFromB.push_back(fromB);
                jPsiMothers.push_back(mothers);
            }
        }

        /// combine D mesons and JPsi mesons
        for (auto ijPsi{0u}; ijPsi < jPsis.size(); ++ijPsi) {
            for (auto iOpenCharm{0u}; iOpenCharm < openCharms.size(); ++iOpenCharm) {
                bool isSPS{false};
                for (auto const& momOpenCharm : openCharmMothers[iOpenCharm]) {
                    if (std::find(jPsiMothers[ijPsi].begin(), jPsiMothers[ijPsi].end(), momOpenCharm) != jPsiMothers[ijPsi].end()) { // we found mother in common
                        isSPS = true;
                        break;
                    }
                }
                tuplePairs->Fill(jPsis[ijPsi].pT(), jPsis[ijPsi].y(), jPsis[ijPsi].eta(), jPsis[ijPsi].phi(), jPsiFromB[ijPsi],
                                 openCharms[iOpenCharm].pT(), openCharms[iOpenCharm].y(), openCharms[iOpenCharm].eta(), openCharms[iOpenCharm].phi(), openCharmFromB[iOpenCharm], openCharms[iOpenCharm].id(), isSPS);
            }
        }
    }
}

//__________________________________________________________________________________________________
void simulateDmesJPsi(int nEvents, int tune, int process, bool usePtHardBins, float energy, int seed, std::string outFileNameRoot)
{
    //__________________________________________________________
    // create and configure pythia generator

    Pythia8::Pythia pythia;
    if(process == kSoftQCD)
    {
        pythia.readString("SoftQCD:inelastic = on"); // we only simulate inelastic processes
    }
    else if(process == kHardQCD)
    {
        pythia.readString("HardQCD:hardccbar = on");
        pythia.readString("HardQCD:hardbbbar = on");
        pythia.readString("HardQCD:gg2ccbar = on");
        pythia.readString("HardQCD:gg2bbbar = on");
        pythia.readString("HardQCD:qqbar2ccbar = on");
        pythia.readString("HardQCD:qqbar2bbbar = on");
    }

    // tune for charmonia
    pythia.readString("CharmoniumShower:all = on");

    // set tune
    if(tune == kMonash)
    {
        pythia.readString(Form("Tune:pp = 14"));
    }
    else if(tune == kCRMode0)
    {
        pythia.readString(Form("Tune:pp = 14"));
        pythia.readString("ColourReconnection:mode = 1");
        pythia.readString("ColourReconnection:allowDoubleJunRem = off");
        pythia.readString("ColourReconnection:m0 = 2.9");
        pythia.readString("ColourReconnection:allowJunctions = on");
        pythia.readString("ColourReconnection:junctionCorrection = 1.43");
        pythia.readString("ColourReconnection:timeDilationMode = 0");
        pythia.readString("StringPT:sigma = 0.335");
        pythia.readString("StringZ:aLund = 0.36");
        pythia.readString("StringZ:bLund = 0.56");
        pythia.readString("StringFlav:probQQtoQ = 0.078");
        pythia.readString("StringFlav:ProbStoUD = 0.2");
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.12");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation = 5");
    }
    else if(tune == kCRMode2)
    {
        pythia.readString(Form("Tune:pp = 14"));
        pythia.readString("ColourReconnection:mode = 1");
        pythia.readString("ColourReconnection:allowDoubleJunRem = off");
        pythia.readString("ColourReconnection:m0 = 0.3");
        pythia.readString("ColourReconnection:allowJunctions = on");
        pythia.readString("ColourReconnection:junctionCorrection = 1.20");
        pythia.readString("ColourReconnection:timeDilationMode = 2");
        pythia.readString("ColourReconnection:timeDilationPar = 0.18");
        pythia.readString("StringPT:sigma = 0.335");
        pythia.readString("StringZ:aLund = 0.36");
        pythia.readString("StringZ:bLund = 0.56");
        pythia.readString("StringFlav:probQQtoQ = 0.078");
        pythia.readString("StringFlav:ProbStoUD = 0.2");
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.15");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation = 5");
    }
    else if(tune == kCRMode3)
    {
        pythia.readString(Form("Tune:pp = 14"));
        pythia.readString("ColourReconnection:mode = 1");
        pythia.readString("ColourReconnection:allowDoubleJunRem = off");
        pythia.readString("ColourReconnection:m0 = 0.3");
        pythia.readString("ColourReconnection:allowJunctions = on");
        pythia.readString("ColourReconnection:junctionCorrection = 1.15");
        pythia.readString("ColourReconnection:timeDilationMode = 3");
        pythia.readString("ColourReconnection:timeDilationPar = 0.073");
        pythia.readString("StringPT:sigma = 0.335");
        pythia.readString("StringZ:aLund = 0.36");
        pythia.readString("StringZ:bLund = 0.56");
        pythia.readString("StringFlav:probQQtoQ = 0.078");
        pythia.readString("StringFlav:ProbStoUD = 0.2");
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.05");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation = 5");
    }

    if (usePtHardBins)
    {
        pythia.readString("PhaseSpace:pTHatMin = 20");
        pythia.readString("PhaseSpace:pTHatMax = 200");
    }

    // init
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed %d", seed));
    pythia.settings.mode("Beams:idA", 2212);
    pythia.settings.mode("Beams:idB", 2212);
    pythia.settings.parm("Beams:eCM", energy);
    pythia.init();

    //__________________________________________________________
    // define outputs
    auto tuplePairs = new TNtuple("tuplePairs", "tuplePairs", "pTJPsi:yJPsi:etaJPsi:phiJPsi:JPsifromB:pTD:yD:etaD:phiD:DfromB:pdgD:commonAnchestor");
    auto histNumEvents = new TH1D("histNumEvents", "", 2, 0., 2.);
    histNumEvents->GetXaxis()->SetBinLabel(1, "num events set");
    histNumEvents->GetXaxis()->SetBinLabel(2, "num events succeded");
    auto histXsec = new TH1D(Form("histXsec_%d", seed), ";cross section (mb)", 1, 0., 1.);
    histXsec->GetXaxis()->SetBinLabel(1, "");
    histNumEvents->SetBinContent(1, nEvents);
    std::map<int, TH3F*> histPtVsY;
    for (auto const& pdg : absPdgOpenCharm) {
        histPtVsY[pdg] = new TH3F(Form("histPtVsY_%d", pdg), ";#it{p}_{T} (GeV/#it{c}); #it{y}; fromB", 500, 0., 50., 200, -10., 10., 2, -0.5, 1.5);
    }
    histPtVsY[443] = new TH3F("histPtVsY_443", ";#it{p}_{T} (GeV/#it{c}); #it{y}; fromB", 500, 0., 50., 200, -10., 10., 2, -0.5, 1.5);
    runSimulation(pythia, tuplePairs, histNumEvents, histPtVsY, nEvents);
    histXsec->SetBinContent(1, pythia.info.sigmaGen());
    histXsec->SetBinError(1, pythia.info.sigmaErr());

    // save root output file
    TFile outFile(outFileNameRoot.data(), "recreate");
    tuplePairs->Write();
    histNumEvents->Write();
    histXsec->Write();
    for (auto const& pdg : absPdgOpenCharm) {
        histPtVsY[pdg]->Write();
    }
    histPtVsY[443]->Write();
    outFile.Close();
}

