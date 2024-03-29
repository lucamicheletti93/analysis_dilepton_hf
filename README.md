# Task configuration for D0 and J/psi tasks

- Run D0 task:
  ```ruby
  o2-analysis-pid-tof-full -b --configuration json://configuration_dq_hf.json | o2-analysis-trackselection -b --configuration json://configuration_dq_hf.json | o2-analysis-multiplicity-table -b --configuration json://configuration_dq_hf.json | o2-analysis-pid-tof-base -b --configuration json://configuration_dq_hf.json | o2-analysis-track-propagation -b --configuration json://configuration_dq_hf.json | o2-analysis-dq-task-jpsi-hf -b --configuration json://configuration_dq_hf.json | o2-analysis-hf-candidate-creator-2prong -b --configuration json://configuration_dq_hf.json | o2-analysis-hf-candidate-selector-d0 -b --configuration json://configuration_dq_hf.json | o2-analysis-timestamp -b --configuration json://configuration_dq_hf.json | o2-analysis-pid-tpc-full -b --configuration json://configuration_dq_hf.json | o2-analysis-track-to-collision-associator -b --configuration json://configuration_dq_hf.json | o2-analysis-pid-tpc-base -b --configuration json://configuration_dq_hf.json | o2-analysis-event-selection -b --configuration json://configuration_dq_hf.json | o2-analysis-hf-track-index-skim-creator -b --configuration json://configuration_dq_hf.json --aod-file AO2D.root
  ```
  
- Run J/psi task:
  ```ruby
  o2-analysis-fwdtrackextension -b --configuration json://configuration_dq_hf.json | o2-analysis-timestamp -b --configuration json://configuration_dq_hf.json | o2-analysis-zdc-converter -b --configuration json://configuration_dq_hf.json | o2-analysis-multiplicity-table -b --configuration json://configuration_dq_hf.json | o2-analysis-event-selection -b --configuration json://configuration_dq_hf.json | o2-analysis-dq-table-maker -b --configuration json://configuration_dq_hf.json | o2-analysis-dq-table-reader -b --configuration json://configuration_dq_hf.json | o2-analysis-track-propagation -b --configuration json://configuration_dq_hf.json --aod-file AO2D.root
  ```

- Run J/psi - D0 task:
  ```ruby
  o2-analysis-pid-tof-full -b --configuration json://configuration_dq_hf.json | o2-analysis-trackselection -b --configuration json://configuration_dq_hf.json | o2-analysis-multiplicity-table -b --configuration json://configuration_dq_hf.json | o2-analysis-pid-tof-base -b --configuration json://configuration_dq_hf.json | o2-analysis-fwdtrackextension -b --configuration json://configuration_dq_hf.json | o2-analysis-timestamp -b --configuration json://configuration_dq_hf.json | o2-analysis-zdc-converter -b --configuration json://configuration_dq_hf.json | o2-analysis-event-selection -b --configuration json://configuration_dq_hf.json | o2-analysis-dq-table-maker -b --configuration json://configuration_dq_hf.json | o2-analysis-dq-table-reader -b --configuration json://configuration_dq_hf.json | o2-analysis-dq-task-jpsi-hf -b --configuration json://configuration_dq_hf.json | o2-analysis-hf-task-d0 -b --configuration json://configuration_dq_hf.json | o2-analysis-hf-candidate-creator-2prong -b --configuration json://configuration_dq_hf.json | o2-analysis-hf-candidate-selector-d0 -b --configuration json://configuration_dq_hf.json | o2-analysis-timestamp -b --configuration json://configuration_dq_hf.json | o2-analysis-pid-tpc-full -b --configuration json://configuration_dq_hf.json | o2-analysis-track-to-collision-associator -b --configuration json://configuration_dq_hf.json | o2-analysis-pid-tpc-base -b --configuration json://configuration_dq_hf.json | o2-analysis-hf-track-index-skim-creator -b --configuration json://configuration_dq_hf.json  | o2-analysis-track-propagation -b --configuration json://configuration_dq_hf.json --aod-file AO2D.root
  ```
  
Test sample:
```ruby
/alice/data/2022/LHC22m/523308/apass4_tpc_v1/0050/o2_ctf_run00523308_orbit0329613340_tf0000082184_epn072/002/AO2D.root
```

Helper script to generate another script (`run_workflow.sh`) that runs the workflows in parallel on multiple inputs
```ruby
chmod +x generate_workflow.sh
./generate_workflow.sh
```
To run the workflows:
```ruby
./run_workflow.sh
```

Workflow to run the task on reduced AO2Ds and produce analysis trees:
```ruby
o2-analysis-dq-task-j-psi-hf -b --configuration json://configuration.json --aod-file @input_data.txt --aod-writer-json OutputDirector.json --aod-parent-access-level 1
```

