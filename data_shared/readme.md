### Links to the shared cernbox directories
- LHC24 pass1 skimmed: https://cernbox.cern.ch/s/AyUkrcuke5ETbfA
- LHC24 pass1 minimum bias (D<sup>0</sup> only): https://cernbox.cern.ch/s/f6lhMNhEjv7zNWR

## MC outputs for D<sup>0</sup>
- LHC24k3 (anchored to LHC24 pass1): https://cernbox.cern.ch/s/crhwnFFlZAv6Ze7

## D0 quantities
- efficiency map: https://github.com/lucamicheletti93/analysis_dilepton_hf/blob/main/data_shared/axeD0.root
- reflection template: https://github.com/lucamicheletti93/analysis_dilepton_hf/blob/main/data_shared/reflections_d0.root
- prompt fraction: https://github.com/lucamicheletti93/analysis_dilepton_hf/blob/main/macro/inclusive_D0/cutvariation/promptfrac_dzero_pp13dot6tev_LHC24_JPsiD_y06.root
- D<sup>0</sup> cross section for p<sub>T</sub> > 0.5 GeV/c and |y|<0.5: https://github.com/lucamicheletti93/analysis_dilepton_hf/blob/main/data_shared/dzero_xsec_pp13dot6TeV.root
- D<sup>0</sup> cross section for p<sub>T</sub> > 0.5 GeV/c and |y|<0.6 : https://github.com/lucamicheletti93/analysis_dilepton_hf/blob/main/data_shared/dzero_xsec_pp13dot6TeV_y06.root

### Useful commands

Inside `fDiMuon` / `fDiElectron` repository run the command:
```ruby
o2-analysis-dq-task-j-psi-hf -b --configuration json://configuration.json --aod-file @input_data.txt --aod-writer-json OutputDirector.json
```
to produce the analysis trees
