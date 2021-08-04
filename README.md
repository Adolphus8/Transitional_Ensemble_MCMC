# Transitional Ensemble Markov Chain Monte Carlo
This repository presents a collection of tutorials (written in MATLAB) which seeks to demonstrate the implementation of the Transitional Ensemble Markov Chain Monte Carlo (TEMCMC) based on the literature by [Lye et. al (2022)](). Currently, 3 tutorials are presented here aimed at allowing users who are new such sampler find their footing around its concept and its workings. The details to these tutorials are as follows:

## Tutorials:

### 1) Illustrative Example:
See: Illustrative_Example.m

This example seeks to provide an illustration of the Affine-invariant property of this Ensemble sampler. 

In most instances, when one encounters a highly-anisotropic target distribution, one approach would be to perform an Affine transformation so as to simplify the form of the target distribution and allow for it to be sampled from easily. From there, an inverse Affine transformation is performed on the generated samples so as to obtain the associated samples that is obtained indirectly from the original highly-anisotropic target distribution. 

For an Affine-invariant sampler, such as the Ensemble sampler here, it is able to sample directly from both the highly-anisotropic and the Affine-transformed distributions and views both distributions as equal. What this means is that an Affine-invariant sampler's sampling performance is unaffected regardless of whether that target distribution is scaled or not. This tutorial aims to highlight this property of the Ensemble sampler and as a comparison, the [Metropolis-Hastings](https://doi.org/10.1093/biomet/57.1.97) sampler will also be implemented to demonstrate the latter's absence of the Affine-invariant property.

### 2) Example Coupled Oscillator:
See: example_CoupledOscillator.m

This tutorial seeks to compare the robustness and strength of the TEMCMC sampler against the standard [TMCMC](https://doi.org/10.1061/(ASCE)0733-9399(2007)133:7(816)) sampler in identifying the spring stiffness of the primary and secondary springs of the Coupled Oscillator as well as the associated noise with the measured frequency modes.

### 3) Example Himmelblau's Function:
See: example_himmelblau.m

This tutorial seeks to compare the effectiveness of the TEMCMC sampler against the standard [TMCMC](https://doi.org/10.1061/(ASCE)0733-9399(2007)133:7(816)) sampler in identifying the peaks of the 4-peaked Himmelbau's function.

## Reference(s):
* A. Lye, A. Cicirello, and E. Patelli (2022). An efficient and robust sampler for Bayesian inference: Transitional Ensemble Markov Chain Monte Carlo. *Mechanical Systems and Signal Processing, 133*(7), 816-832. doi: 

## Author:
* Name: Adolphus Lye
* Contact: adolphus.lye@liverpool.ac.uk
* Affiliation: Insitute for Risk and Uncertainty, University of Liverpool
