%FIGURE_MORE_INITS Showcase the (absence of) effect of initial conditions

clear; clc; close all;
init_prj;


%% Get motifs with 120 random initial conditions

mdlName = 'HCP_Rest_FIX_Simple_Mdl200_sess';
outFile = fullfile('data', ['motifs_randInit_' mdlName '.mat']);
nSubs = 0;

% Get motifs
Report_Attractor_Landscape(mdlName, 'outFile', outFile, 'vis', false, ...
    'nSubs', nSubs, 'init_w_dat', false, 'nInits', 120);

% Load the motifs
load(outFile, 'allEq', 'allLCMotif');
allEq2 = allEq;
allLCMotif2 = allLCMotif;
clear allEq allLCMotif;
if nSubs == 0
    nSubs = size(allEq2, 1);
end


%% Load the motifs obtained using 1000 samples from deconvolved data

load(fullfile('data', 'motifs_HCP_Rest_FIX_Simple_Mdl200_sess.mat'), ...
    'allEq', 'allLCMotif', 'censored');
allEq = allEq(1:nSubs, :);
allLCMotif = allLCMotif(1:nSubs, :);
censored = censored(1:nSubs, :);
censored = any(censored, 2);

% Remove the censored models
allEq(censored, :) = [];
allLCMotif(censored, :) = [];
allEq2(censored, :) = [];
allLCMotif2(censored, :) = [];


%% Compare the motifs

% Number and type
nEq = cellfun(@(x) size(x, 2), allEq);
nEq2 = cellfun(@(x) size(x, 2), allEq2);
nLC = cellfun(@(x) size(x, 2), allLCMotif);
nLC2 = cellfun(@(x) size(x, 2), allLCMotif2);
idx = (nEq == nEq2 & nLC == nLC2);
fprintf('Proportion of models with number and type of motifs matching: %.2f%%\n', ...
    mean(idx, "all") * 100);

% Focus on the models with matching number and type of motifs
eq1 = allEq(idx(:));
eq2 = allEq2(idx(:));
lc1 = allLCMotif(idx(:));
lc2 = allLCMotif2(idx(:));

% Pick the first motif from every limit cycle
for i = 1:numel(lc1)
    if ~isempty(lc1{i})
        % Add FlipSign to enforce the direction
        tmp = cellfun(@(x) FlipSign(x(:, 1)), lc1{i}, 'UniformOutput', false);
        lc1{i} = cell2mat(tmp);
        tmp = cellfun(@(x) FlipSign(x(:, 1)), lc2{i}, 'UniformOutput', false);
        lc2{i} = cell2mat(tmp);
    else
        lc1{i} = [];
        lc2{i} = [];
    end
end

% Calculate the distance between the motifs
eqDist = cellfun(@(x, y) pdist2(x', y'), eq1, eq2, 'UniformOutput', false);
lcDist = cellfun(@(x, y) pdist2(x', y'), lc1, lc2, 'UniformOutput', false);

% Find optimal matching mean distance
totalDist = zeros(size(eqDist));
for i = 1:numel(eqDist)
    if ~isempty(eqDist{i})
        [~, cost] = assignmentoptimal(eqDist{i});
        totalDist(i) = totalDist(i) + cost;
    end
    if ~isempty(lcDist{i})
        [~, cost] = assignmentoptimal(lcDist{i});
        totalDist(i) = totalDist(i) + cost;
    end
    totalDist(i) = totalDist(i) / (size(eqDist{i}, 1) + size(lcDist{i}, 1));
end

% Normalized by mean Euclidean norm of motifs
motifNorm = zeros(size(eqDist));
for i = 1:numel(eqDist)
    if ~isempty(eqDist{i})
        motifNorm(i) = motifNorm(i) + sum(vecnorm(eq1{i}));
    end
    if ~isempty(lcDist{i})
        motifNorm(i) = motifNorm(i) + sum(vecnorm(lc1{i}));
    end
    motifNorm(i) = motifNorm(i) / (size(eqDist{i}, 1) + size(lcDist{i}, 1));
end


%% Plot the results

figure('Units', 'inches', 'Position', [1 1 9 6]);
tlo = tiledlayout(1, 2, "Padding", "compact");

tmp = [mean(idx(:)), 1 - mean(idx(:))] * 100;
nexttile;
hold on;
bar(tmp);
ylabel('Proportion (%)');
xticks(1:2);
xticklabels({'Matched', 'Mismatched'});
text(1:2, tmp, compose('%.1f', tmp), "HorizontalAlignment", "center", ...
    "VerticalAlignment", "bottom");
ylim([0 105])
fontsize(gca, 12, "points");
title('Number and type of attractors', 'FontSize', 14);
subtitle(sprintf('n = %d', numel(allEq)));

nexttile;
histogram(log10(totalDist ./ motifNorm), -16:0);
xticks(-16:2:0)
xlabel('Log10(mean distance / mean norm)');
ylabel('Frequency');
fontsize(gca, 12, "points");
title({'Normalized distance between', 'matched motifs'}, 'FontSize', 14);

% title(tlo, 'Comparing motifs extracted using different initial conditions', ...
%     'FontSize', 16);
PrintAsSeen(gcf, fullfile('figures', mdlName, 'more_inits'), '-dpng', '-r300')