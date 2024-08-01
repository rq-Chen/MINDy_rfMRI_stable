%FIGURE_ROBUST_KMEANS Use robust K-means to repeat the clustering

init_prj;
clear; clc; close all;


%% Parameters

mdlName = 'HCP_Rest_FIX_Simple_Mdl200_sess';
mtf_csv = fullfile('data', ['motifs_' mdlName, '.csv']);
kmeans_csv = fullfile('data', ['kmeans_' mdlName, '.csv']);


%% Load data and save to csv

% Get the vectors to cluster
[V, VNames, clIdx, C] = Report_Motifs(mdlName, 'vis', false, 'Klist', 4, ...
    'normalize', true);  % Scale each motif to unit length

% Convert V from (nSubs, nSess) cell of (nParcels, nij) matrix to
% (nMotifs, nParcels) matrix
V = [V{:}]';  % (nMotifs, nParcels) matrix
VNames = [VNames{:}]';  % (nMotifs, 1) string
clIdx = [clIdx{:}]';  % (nMotifs, 1) int

% Get parcel names
atlas = GetYeoNetworks(size(V, 2));

% Save to csv
writetable(array2table(V, 'VariableNames', atlas.ParcelNames), mtf_csv);
writetable(array2table(clIdx, 'VariableNames', {'Cluster'}), kmeans_csv);


%% Run robust K-means in R

% system('Rscript code/Supplementary/Figure_Robust_Kmeans.R');


%% Load results and plot

tbl = readtable(fullfile('data', ['robust_kmeans_k4_centers_' mdlName '.csv']));
C = tbl{2:end, 2:end};
tmp = max(abs(C(:)));
lim = [-tmp tmp];
f = figure('Units', 'inches', 'Position', [1 1 6 6]);
tlo = tiledlayout(2, 2, "TileSpacing", "tight", "Padding", "tight");
for i = 1:4
    nexttile;
    PlotYeoSurface(C(i, :));
    clim(lim);
    title(sprintf('Cluster %d', i), 'FontSize', 14)
end
cb = colorbar;
cb.Layout.Tile = 'east';
PrintAsSeen(fullfile('figures', mdlName, ['Robust_Kmeans_K4_' mdlName]), '-dpng', '-r300');