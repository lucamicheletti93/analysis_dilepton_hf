### Links to the shared cernbox directories
- LHC23 pass4 skimmed: https://cernbox.cern.ch/s/DYW65twa2sWH8qz

### Useful commands

Inside `fDiMuon` / `fDiElectron` repository run the command:
```ruby
o2-analysis-dq-task-j-psi-hf -b --configuration json://configuration.json --aod-file @input_data.txt --aod-writer-json OutputDirector.json
```
to produce the analysis trees
