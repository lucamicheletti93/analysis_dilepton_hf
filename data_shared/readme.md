### Links to the shared cernbox directories
- LHC22 pass7 skimmed: https://cernbox.cern.ch/s/nJ8YNoR77FIDwvJ
- LHC23 pass4 skimmed: https://cernbox.cern.ch/s/DYW65twa2sWH8qz
- LHC24 pass1 skimmed: https://cernbox.cern.ch/s/gdXXO4gc1YBKlNX

## MC outputs for D<sup>0</sup>
- LHC24g5 (anchored to LHC22 pass7): https://cernbox.cern.ch/s/n6fglJh0XOpbWni
- LHC24h1 (anchored to LHC23 pass4): https://cernbox.cern.ch/s/HbZVGqFlnbHrQTU

### Useful commands

Inside `fDiMuon` / `fDiElectron` repository run the command:
```ruby
o2-analysis-dq-task-j-psi-hf -b --configuration json://configuration.json --aod-file @input_data.txt --aod-writer-json OutputDirector.json
```
to produce the analysis trees
