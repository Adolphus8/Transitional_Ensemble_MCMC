# Transitional Ensemble Markov Chain Monte Carlo
This repository presents a collection of tutorials (written in MATLAB) which seeks to demonstrate the implementation of the Transitional Ensemble Markov Chain Monte Carlo (TEMCMC) based on the literature by [Lye et. al (2022)](https://doi.org/10.1016/j.ymssp.2021.108471). Currently, 3 tutorials are presented here aimed at allowing users who are new such sampler find their footing around its concept and its workings. The details to these tutorials are as follows:

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

### 4) Aluminium Frame Problem:
See: Aluminium_Frame_Problem folder

This problem involves a 2 Degree-of-Freedom Shear Alumium Frame with 2 moveable masses, whose respective positions represent the location of damage. Here, actual experimental data of the response frequencies of the frame is obtained from a hammer test given a configuration of the positions of the moveable masses. This data is then used to infer the mass positions using an Artificial Neural Network (ANN), used in [R. Rocchetta et. al (2018)](https://doi.org/10.1016/j.ymssp.2017.10.015), that is trained with 200 synthetic data of the output response frequencies and the input mass positions obtained from the frame's Finite Element Model and Monte Carlo sampling.

For this problem, the TEMCMC and TMCMC samplers are used to perform robust Bayes with different aggregated likelihood functions and obtain Posteriors of the mass positions in the form of P-boxes from which the interval estimates of the mass positions are obtained. The purpose of this work is to verify the estimates by the TEMCMC and validate the algorithm using actual experimental data.

For more details to the Alumium Frame set-up, readers can also refer to the work by [H. H. Khodaparast et. al (2011)](https://doi.org/10.1016/j.ymssp.2010.10.009), [P. Liang et. al (2016)](http://past.isma-isaac.be/downloads/isma2016/papers/isma2016_0351.pdf), and [Z. Yuan et. al (2019)](https://doi.org/10.1016/j.ymssp.2018.05.048).

## Reference(s):
* A. Lye, A. Cicirello, and E. Patelli (2022). An efficient and robust sampler for Bayesian inference: Transitional Ensemble Markov Chain Monte Carlo. *Mechanical Systems and Signal Processing, 167*, 108471. doi: [10.1016/j.ymssp.2021.108471](https://doi.org/10.1016/j.ymssp.2021.108471)
* A. Lye, A. Cicirello, and E. Patelli (2021). Sampling methods for solving Bayesian model updating problems: A tutorial. *Mechanical Systems and Signal Processing, 159*, 107760. doi: [10.1016/j.ymssp.2021.107760](https://doi.org/10.1016/j.ymssp.2021.107760)

* H. H. Khodaparast, J. E. Mottershead, and K. J. Badcock (2011). Interval model updating with irreducible uncertainty using the Kriging predictor. *Mechanical Systems and Signal Processing, 25*(4), 1204-1226. doi: [10.1016/j.ymssp.2010.10.009](https://doi.org/10.1016/j.ymssp.2010.10.009)
* P. Liang, J. E. Mottershead, and F. A. DiazDelaO (2016). Model Updating with the Kriging Predictor: Effect of Code Uncertainty. [*In the Proceedings of ISMA 2016 including USD 2016*](http://past.isma-isaac.be/downloads/isma2016/papers/isma2016_0351.pdf)
* R. Rocchetta, M. Broggi, Q. Huchet, and E. Patelli (2018). On-line Bayesian model updating for structural health monitoring. *Mechanical Systems and Signal Processing, 103*, 174-195. doi: [10.1016/j.ymssp.2017.10.015](https://doi.org/10.1016/j.ymssp.2017.10.015)
* Z. Yuan, P. Liang, T. Silva, K. Yu, and J. E. Mottershead (2019). Parameter selection for model updating with global sensitivity analysis. *Mechanical Systems and Signal Processing, 115*, 483-496. doi: [10.1016/j.ymssp.2018.05.048](https://doi.org/10.1016/j.ymssp.2018.05.048)

## Author:
* Name: Adolphus Lye
* Contact: adolphus.lye@liverpool.ac.uk
* Affiliation: Insitute for Risk and Uncertainty, University of Liverpool
