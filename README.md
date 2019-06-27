# Markov Melding
ST955: Dissertation, MSc Statistics, University of Warwick

## Meeting 26/06/2019
To-do:
- [ ] Optimise sampler via acceptance ratio or ESJD on one of the components then use covariance of `mle` to find the scaling for other components. *Edit: not sure that this is a good idea, or maybe my optimisation procedure is just a bit janky*
- [x] Facilitate above by adding acceptance ratio output to `mwg` function potentially
- [ ] "Meld" models 1 and 2 (with consistent priors first if this makes sense, then with inconsistent priors)