# MINDy Analysis

This repo contains the analysis scripts for:

Chen, R., Singh, M., Braver, T.S., & Ching, S. (2025). Dynamical models reveal anatomically reliable attractor landscapes embedded in resting state brain networks. Imaging Neuroscience, Advance Publication. [https://doi.org/10.1162/imag_a_00442](https://doi.org/10.1162/imag_a_00442)

For simplicity, we include some functions from other packages in this repo. Please note that they might be covered by different licenses. In particular, we include `MINDy_Base_v1.0` from [singhmf/MINDy](https://github.com/singhmf/MINDy) (and modified it a little bit), `ft_cifti` (`cifti-matlab-master`) from [Washington-University/cifti-matlab](https://github.com/Washington-University/cifti-matlab), and part of `dFCwalk` from [FunDyn/dFCwalk](https://github.com/FunDyn/dFCwalk).

The scripts were tested using MATLAB R2022b under Windows and Linux. One Supplementary Figure requires `R` and the package `wrsk` for robust K-means clustering.

## Pipeline

1. Follow [Singh2020PreProc](https://github.com/rq-Chen/Singh2020PreProc) to preprocess the HCP rfMRI data. By default there will be one file for each subject under the output folder. Files are named like "sub100206Y02.mat" (200 Schaefer atlas parcels + 19 freesurfer subcortical ROIs, for participant 100206).

2. Run `init_prj.m` to set up paths (see below).

3. Run `code/training/GetHCPRestModel.m` to train the MINDy models and save the models as well as the deconvoluted rfMRI data into `data/`.

4. Run `code/Report_Attractor_Landscape.m` to visualize the vector field of every single model and obtain the attractors. Figures (128 in total!) will be saved to `figures/` and attractors will be saved to `data/`.

5. Visually inspect the figures to make sure the attractors were correctly identified. Record the models where the numerical procedure might have failed in the file `data/censored_models.csv`.

6. Run `code/Censor_Models.m` to exclude the problematic models (as well as the models from the same participant) from further analysis and re-generate the figure of the distribution of dynamics across clean models.

7. Run `code/Report_Bifurcation.m` to perform the induced bifurcation analysis.

8. Run `code/Report_Attractor_Reliability.m` to perform the attractor reliability analysis.

9. Run `code/Report_Motifs.m` to perform the attractor clustering analysis.

10. Run other scripts under `code/Supplementary/` to generate the supplementary figures and videos (see below).

## Figures

- 1B: From `code/Report_Motifs.m`.

- 1C: From `code/example_dynamics.m` with `plotType = 'popular'`.

- 2: From `code/Report_Bifurcation.m`.

- 3: From `code/Report_Attractor_Reliability.m`.

- 4: From `code/Report_Motifs.m`.

- S1: `code/Supplementary/Figure_Cross_Validation.m`.

- S2: `code/Supplementary/Report_Simulation_Validation.m`.

- S3: `code/Supplementary/Figure_Shuffle_Weights.m` and `code/Supplementary/Figure_Shuffle_Weights.Rmd`.

- S4: `code/Supplementary/Report_Analysis_Deconv.m`

- S5: `code/Supplementary/Report_Param_Reliability.m`.

- S6: `code/Supplementary/MINDy200_CCA.m` (Note: CCA was performed by [this repo](https://github.com/rq-Chen/HCP_CCA_1200_OSF)).

- S7: `code/Supplementary/Figure_SC.m`

- S8: `code/Censor_Models.m`

- S9: `code/Supplementary/Figure_more_inits.m`

- S10: `code/Censor_Models.m`

- S11: `code/example_dynamics.m` with `plotType = 'rare'`.

- S12: `code/Supplementary/infinite_period_bifurcation.m`.

- S13: `code/Report_Motifs.m`.

- S14-S16: `code/Report_Motifs.m` with corresponding inputs.

- S17: `code/Supplementary/Figure_Robust_Kmeans.m` and `code/Supplementary/Figure_Robust_Kmeans.R`

- S18-19: `code/Report_Motifs.m`.

- S20-22: Similar to S10, S13 and Figure 4, but using models fit with one-step instead of two-step difference. This is done by setting `'MINDyType, 'NoSmooth'` when running `GetHCPRestModel`.

- Videos: `code/Supplementary/Video_trajectories.m`.

## Directory and file structure

### Paths

In most scripts, we encode the paths relative to the project root, so it's required to set the working directory as project root (where the git folder is in). Besides, we need to add the scripts to MATLAB search path. You can run the script `init_prj.m` to set up `pwd` and the search path.

### Data

We saved the deconvoluted rfMRI data to the `data` folder. Using default parameters, the file will be called `HCP_Rest_FIX_Simple_Deconv200_sess.mat`. The file contains:

- an `(nSubj, nRuns)` cell array `allDat`: the deconvoluted rfMRI data. Each cell is a `(nParcels, nTRs)` matrix (TR = .72 second). The data is preprocessed, parcellated, deconvoluted, smoothed and z-scored (ready for MINDy fitting or comparison with simulations).

- an `(nSubj, nModels)` cell array `runIdx`: the indices of the runs used to train each model. By default, `nModels = 2` (one model for each session) and `runIdx(:, 1) = {1:2}` and `runIdx(:, 2) = {3:4}`.

- some other variables like subject IDs, input files.

### Models

The MINDy models were also saved to the `data` folder. Using default parameters, the file will be called `HCP_Rest_FIX_Simple_Mdl200_sess.mat`. The file contains similar variables as the deconvoluted rfMRI data, but with the models instead of the data. The models are saved as a `(nSubj, nModels)` cell array `allMdl`. Each cell is a structure. In particular, the `Param` field contains the parameters of the model. It is a `(1, 6)` cell:

- $W_s$, the sparse component of the connectivity matrix, `(nParcels, nParcels)`
- $\alpha$, curvature of the transfer function of the neural mass model, `(nParcels, 1)`
- a 1 \* 2 matrix:
  - $b$, slope of the transfer function, currently fixed as $20/3$
  - a shift term for the transfer function, currently fixed as $0$
- $c$, a constant drive term, curretly fixed as $0$
- $W$, connectivity matrix, sum of $W_s$ (sparse component) and $W_l$ (low-rank component), `(nParcels, nParcels)`
- $D$, decay coefficient, `(nParcels, 1)`

For most uses, `W = Param{5}; D = Param{6}; A = Param{2}`.

### Attractors

The attractors were also saved to the `data` folder. Using default parameters, the file is called `motifs_PC_HCP_Rest_FIX_Simple_Mdl200_sess.mat`. It contains:

- `allEq`: `(nSubj, nModels)` cell array of `(nParcles, ?)` matrices where each column is a stable equilibrium (`?` could be zero).
- `allLCMotifs`: `(nSubj, nModels)` cell array of `(1, ?)` cell array (one cell for each limit cycle, `?` could be zero). Each cell contains one `(nParcels, 4)` matrix, representing the "positive/negative extremes" of the "major/minor axes" of the limit cycle. The exact definition depends on `code/utilities/LC2Motif.m`. By default, the first column will be the slowest point on the limit cycle, and the rest will be the points closest to 90/180/270 degree to the first point after projecting the limit cycle into first two PCs. Only the first point (i.e., slowest point) is used in our analysis.