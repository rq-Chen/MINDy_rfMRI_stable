%%EXAMPLE_DYNAMICS - Example vector fields for popular and rare dynamics

clear; clc; close all;
init_prj;


%% Parameters

plotType = 'rare';
[types, subjID, iSess] = GetPlotType(plotType);


%% Load data

load(fullfile('data', 'HCP_Rest_FIX_Simple_Mdl200_sess.mat'), 'allMdl', 'sublist');
sublist = extractBetween(sublist, 'sub', 'Y02.mat');
iSub = arrayfun(@(x) find(strcmp(sublist, x), 1), subjID);


%% Plot

N = numel(iSub);
f = figure('Units', 'inches', 'Position', [0 0.5 6 4]);
t = tiledlayout(2, 2, 'TileSpacing', 'compact');
for i = 1:N
    mdl = allMdl{iSub(i), iSess(i)};
    [eq, lc, sims] = GetEqLC(mdl);
    mtf = LC2Motif(lc);
    nexttile();
    lgd = PlotAttractorLandscape(sims(:, 15:5:end, 1:40), eq, mtf, 3, 0);
    title(types{i});
end
lgd.Layout.Tile = 'east';
lgd.Visible = 'on';
lgd.FontName = 'helvetica';
title(t, 'Example vector fields');
PrintAsSeen(f, fullfile('figures', 'HCP_Rest_FIX_Simple_Mdl200_sess', ...
    [plotType '_dynamics']), '-dpdf', '-vector');
PrintAsSeen(f, fullfile('figures', 'HCP_Rest_FIX_Simple_Mdl200_sess', ...
    [plotType '_dynamics']), '-dsvg', '-vector');