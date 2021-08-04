# Affine Invariant Ensemble Sampler Tutorials
This repository presents a collection of tutorials (written in MATLAB) which seeks to demonstrate the implementation of the Affine Invariant Ensemble Sampler based on the works by [Goodman and Weare (2010)](https://projecteuclid.org/euclid.camcos/1513731992). The associated literature is available in the "*References*" folder. Currently, 3 tutorials are presented here aimed at allowing users who are new such sampler find their footing around its concept and its workings. To help facilitate the users in understanding the workings of the sampler as well as the objectives of the tutorials in a simplified manner, a summarised set of notes is also made availale in the "*References*" folder. The details to these tutorials are as follows:

## Tutorials:

### 1) Example Comparison:
This tutorial seeks to provide an illustration of the Affine-invariant property of this Ensemble sampler. 

In most instances, when one encounters a highly-anisotropic target distribution, one approach would be to perform an Affine transformation so as to simplify the form of the target distribution and allow for it to be sampled from easily. From there, an inverse Affine transformation is performed on the generated samples so as to obtain the associated samples that is obtained indirectly from the original highly-anisotropic target distribution. 

For an Affine-invariant sampler, such as the Ensemble sampler here, it is able to sample directly from both the highly-anisotropic and the Affine-transformed distributions and views both distributions as equal. What this means is that an Affine-invariant sampler's sampling performance is unaffected regardless of whether that target distribution is scaled or not. This tutorial aims to highlight this property of the Ensemble sampler and as a comparison, the [Metropolis-Hastings](https://doi.org/10.1093/biomet/57.1.97) sampler will also be implemented to demonstrate the latter's absence of the Affine-invariant property.

### 2) Example Rosenbrock:
This tutorial seeks to evaluate the robustness and strength of the Ensemble sampler in sampling from a highly-anisotropic distribution which cannot be re-scaled via Affine-transformation. As a comparison, the [Metropolis-Hastings](https://doi.org/10.1093/biomet/57.1.97) sampler will also be implemented to demonstrate the relative strength of the Ensemble sampler over the latter in samplign from such complex distributions.

### 3) Example Bi-modal Posterior:
This tutorial seeks to evaluate the performance of the Ensemble sampler in sampling from a bi-modal target distribution. As a comparison, the established [Transitional Markov Chain Monte Carlo (TMCMC)](https://doi.org/10.1061/(ASCE)0733-9399(2007)133:7(816)) sampler will also be implemented as a benchmark/comparison. 

## References:
* J. Goodman and J. Weare (2010). Ensemble Samplers with Affne Invariance.
*Communications in Applied Mathematics and Computational Science, 5*(1), 65-80. doi: 10.2140/camcos.2010.5.65

* W. K. Hastings, W. K. (1970). Monte Carlo sampling methods using Markov chains and their applications. *Biometrika, 57*(1), 97-109. doi: 10.1093/biomet/57.1.97

* H. H. Rosenbrock (1960). An Automatic Method for Finding the Greatest or Least Value of a Function. *The Computer Journal, 3*(3), 175-184. doi: 10.1093/comjnl/3.3.175

* J. Ching, and Y. Chen (2007). Transitional Markov Chain Monte Carlo Method for Bayesian Model Updating, Model Class Selection, and Model Averaging. *Journal of Engineering Mechanics, 133*(7), 816-832. doi: 10.1061/(asce)0733-9399(2007)133:7(816) 

## Author:
* Name: Adolphus Lye
* Contact: adolphus.lye@liverpool.ac.uk
* Affiliation: Insitute for Risk and Uncertainty, University of Liverpool

## Acknowledgement:
The Ensemble sampler code presented here is adopted and modified from the original version (in MATLAB as well) written by Aslak Grinsted. The repository of his work can be accessed [HERE](https://github.com/grinsted/gwmcmc). 
