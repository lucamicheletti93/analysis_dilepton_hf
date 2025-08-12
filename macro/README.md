# Fit code: how to use

Run J/&#968; - D0 fit with the command:
  ```ruby
  python fit.py configs/config_fit_fwd.yml --do_fit
  ```
other options:
- If you want to apply selections to the input sample run:
  ```ruby
  python fit.py configs/config_fit_fwd.yml --do_prefilter
  ```
  A new sample with the given suffix will be created
- If you want to fit the J/&#968; and D0 distributions to extract their shape run:
  ```ruby
  python fit.py configs/config_fit_fwd.yml --do_prefit
  ```
- If you want to apply weights according to a 2D Axe map you have to run:
  ```ruby
  python fit.py configs/config_fit_fwd.yml --do_weightdata
  ```
  A new sample with the suffix `_weighted` will be created. Pay attention to the name of the output tree

# Multi-trial fit: how to use

- Run J/&#968; - D0 multi-trial fit with the command:
  ```ruby
  python multi_trial.py configs/multi_trial/config_fit_CB2_Cheby_CB_Cheby_data_tails.yml --run
  ```
- Combine the systematics using:
  ```ruby
  python multi_trial.py --combine -t multi_trial_list_weighted_fit.txt
  ```
