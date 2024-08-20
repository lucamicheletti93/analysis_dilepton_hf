# Fit code: how to use

Run J/&#968; - D0 fit with the command:
  ```ruby
  python fit.py config_fit.yml --do_fit
  ```
other options:
- If you want to apply selections to the input sample run:
  ```ruby
  python fit.py config_fit.yml --do_prefilter
  ```
  A new sample with the given `suffix` will be created
- If you want to fit the J/$\psi$ and D0 distributions to extract their shape run:
  ```ruby
  python fit.py config_fit.yml --do_prefit
  ```
