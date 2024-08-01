%VIDEO_TRAJECTORIES - Generate some video of simulated trajectories

clear; clc; close all;
init_prj;


%% Parameters

[types1, subjID1, iSess1] = get_example_model_idx('popular');
[types2, subjID2, iSess2] = get_example_model_idx('rare');
types = [types1, types2];
subjID = [subjID1, subjID2];
iSess = [iSess1, iSess2];
clear types1 types2 subjID1 subjID2

mdlName = 'HCP_Rest_FIX_Simple_Mdl200_sess';
figdir = fullfile('figures', mdlName, 'Noisy_Simulation_Gifs');


%% Load data

load(fullfile('data', [mdlName '.mat']), 'allMdl', 'sublist');
sublist = extractBetween(sublist, 'sub', 'Y02.mat');
iSub = arrayfun(@(x) find(strcmp(sublist, x), 1), subjID);


%% Plots

currType = '2FP 0LC';
idx = find(strcmp(types, currType));
dW = 0;
t = 1:10:700;
mdl = allMdl{iSub(idx), iSess(idx)};
fname = fullfile(figdir, ['PCSurf_' sprintf('dW_%02d_', dW * 100) ...
    replace(currType, ' ', '_') '.mp4']);
PlotPCSurfVideo(mdl, fname, dW, t, false);