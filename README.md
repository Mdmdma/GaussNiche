# GaussNiche
A lightweight R framework for designing virtual species in multivariate environmental space using bivariate Gaussian functions fitted on PCA axes, with back-projection to geographic space and built-in diagnostics for pseudo-absence quality (class overlap and sampling bias).

# Methodology 
## 1. Environmental PCA. 
A set of bioclimatic rasters is reduced to its principal components using a correlation-matrix PCA (variables centred and scaled). All subsequent steps operate on the scores of the first two principal components (PC1, PC2), which define the bivariate environmental space available to species within the study area.

## 2. Background characterisation. 
The full set of non-NA raster cells is extracted and used as the background — the universe of environmentally available conditions. Its distribution in PC space is visualised as a kernel density estimate, allowing deliberate placement of the virtual species niche relative to common and rare environments.

## 3. Niche definition. 
The species ecological niche is defined as a bivariate Gaussian function centred on a user-specified optimum (μ) in PC space, with a covariance matrix (Σ) controlling niche breadth along each axis and the correlation between axes. Suitability values are normalised by the theoretical peak of the Gaussian at μ, so the optimum cell always receives a suitability of exactly 1.0, independently of how much background falls near μ.

## 4. Presence sampling.
Each background cell is treated as an independent Bernoulli trial with success probability equal to its suitability value. The resulting binary vector defines the virtual species presence/absence across the landscape, with geographic prevalence emerging naturally from the interaction between niche position and the distribution of available environments.

## 5. Pseudo-absence sampling. 
A set of pseudo-absences is drawn from the background according to a user-specified sampling strategy. The framework is designed to accommodate multiple strategies (random, buffer-based, environmentally uniform) through a modular sampler interface, allowing controlled comparisons between approaches.

## 6. Class overlap. 
The degree of overlap between presences and pseudo-absences in PC space is quantified as the intersection volume of two hypervolumes fitted separately on the two point sets. High overlap indicates that pseudo-absences are drawn from environments similar to those occupied by the species, which can compromise model discrimination.

## 7. Sampling bias. 
The environmental coverage of pseudo-absences is assessed as the ratio between the range of PC1 and PC2 values spanned by the pseudo-absence set and the corresponding background range. Values close to 1 indicate that pseudo-absences span the full environmental gradient available; values substantially below 1 indicate geographic or environmental clustering.

## 8. Geographic back-projection. 
Suitability values, presence/absence, and pseudo-absence locations are mapped back onto the original geographic space by reconstructing rasters from the cell coordinates retained throughout the pipeline.

# Repository Contents

* 1_developing_framework.R: step-by-step development of the full pipeline, from environmental PCA and niche definition to pseudo-absence sampling and diagnostic visualisation. Intended as a transparent illustration of the methodology before abstraction into functions.

* virtualSpecies_fn.R: modular function wrappers encapsulating the core pipeline (niche definition, Bernoulli sampling, pseudo-absence sampling, class overlap, sampling bias diagnostics, geographic back-projection). The pseudo-absence sampler is exposed as a swappable argument to facilitate future comparisons between sampling strategies.

* 2_testing_wrapper_function.R: testing script for the wrapper function, verifying outputs and diagnostic plots across different niche parameter configurations.