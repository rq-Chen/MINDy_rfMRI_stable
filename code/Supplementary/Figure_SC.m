%FIGURE_SC Compare MINDy connectivity matrix and structural connectivity matrix

clear; clc; close all;
init_prj;


%% Load data

mdlName = 'HCP_Rest_FIX_Simple_Mdl200_sess';

% Population-averaged MINDy connectivity matrix W
load(fullfile('data', [mdlName, '.mat']), 'allMdl');
W = cellfun(@(x) x.Param{5}, allMdl, "UniformOutput", false);
W = cat(3, W{:});
W = mean(W, 3);

% Structural connectivity matrix
[SC, ~, tSC, pSC] = GetSCMask();

% Mean functional connectivity matrix
load(fullfile('data', [replace(mdlName, 'Mdl', 'Deconv'), '.mat']), 'allDat', 'runIdx');
FC = cell(size(runIdx));
for i = 1:size(FC, 1)
    for j = 1:size(FC, 2)
        FC{i, j} = corr([allDat{i, runIdx{i, j}}]');
        % Fisher's z-transform
        FC{i, j} = atanh(FC{i, j});
    end
end
FC = mean(cat(3, FC{:}), 3);

% Save the data
save(fullfile('data', 'Population_Averaged_Connectivity.mat'), 'W', 'SC', 'FC', 'tSC', 'pSC');


%% Plot

mat2vec = @(m, varargin) m(triu(true(size(m)), varargin{:}));

Wv = mat2vec(W, 1);
SCv = mat2vec(SC, 1);
FCv = mat2vec(FC, 1);
pSCv = mat2vec(pSC, 1);
mask = pSCv < 0.05 / numel(pSCv);
colors = lines(2);
% c = colors(mask + 1, :);
c = colors(1, :);

f = figure('Units', 'inches', 'Position', [0 0 18 12]);
tlo = tiledlayout(2, 3, 'TileSpacing', 'loose', 'Padding', 'loose');

nexttile();
MyNetMatYeo(FC);
axis square
title('Functional Connectivity Matrix FC', 'FontSize', 14);

nexttile();
MyNetMatYeo(W, [], [], 'none');
axis square
title('MINDy Connectivity Matrix W', 'FontSize', 14);

nexttile();
MyNetMatYeo(SC, [], [], 'none');
axis square
title('Structural Connectivity Matrix SC', 'FontSize', 14);

nexttile();
ScatterSmooth(Wv, FCv, 1, c, 'filled');
xlabel('W');
ylabel('FC');
fontsize(gca, 14, 'points');
title('Correlation between W and FC', 'FontSize', 14);

nexttile();
ScatterSmooth(log(SCv), FCv, 1, c, 'filled');
xlabel('log(SC)');
ylabel('FC');
fontsize(gca, 14, 'points');
title('Correlation between log(SC) and FC', 'FontSize', 14);

nexttile();
ScatterSmooth(log(SCv), Wv, 1, c, 'filled');
xlabel('log(SC)');
ylabel('W');
fontsize(gca, 14, 'points');
title('Correlation between log(SC) and W', 'FontSize', 14);

PrintAsSeen(f, fullfile('figures', mdlName, 'Connectivity_Comparison'), '-dpdf', '-vector');