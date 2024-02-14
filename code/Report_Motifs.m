%% Report_Motifs.m - Analysis of motif clusters
function [V, VNames, clIdx, C] = Report_Motifs(mdlName, varargin)
%REPORT_MOTIFS Analysis of motif clusters
%
%   Usage:
%       [V, VNames, clIdx, C] = Report_Motifs(mdlName, varargin)
%
%   Positional input:
%       mdlName - Name of the model (default: 'HCP_Rest_FIX_Simple_Mdl200_sess')
%
%   Key-value inputs:
%       'mtfFile' - Path to the motif file (default: 'data/motifs_<mdlName>.mat')
%       'vis' - Whether to visualize the results (default: true)
%       'figdir' - Path to the figure directory (default: 'figures/<mdlName>')
%       'mtfType' - Type of motifs to cluster (default: 'all')
%           'FP' - Fixed points
%           'LC' - Limit cycles
%           'all' - All motifs
%       'exclude_origin' - Whether to exclude the stable origin (default: true)
%       'normalize' - Whether to normalize the length of motifs (default: true)
%       'demean' - Whether to subtract the global mean from each parcel (default: false)
%       'Klist' - List of numbers of clusters to try (default: 2:10)
%       'K' - Number of clusters to use (default: 4)
%
%   Outputs:
%       V - (nSubs, nSess) cell array of (nParcels, nij) motifs
%       VNames - (nSubs, nSess) cell array of (1, nij) strings, motif names
%       clIdx - (nSubs, nSess) cell array of (1, nij) cluster indices
%       C - (K, nParcels) cluster centroids
%
%   Note that 'normalize' will not actually influence the clustering except
%   scaling the results, since we use cosine distance for clustering.


%% Parameters

if nargin == 0 || isempty(mdlName)
    mdlName = 'HCP_Rest_FIX_Simple_Mdl200_sess';
end

p = inputParser;
addParameter(p, 'mtfFile', fullfile('data', ['motifs_' mdlName '.mat']), @ischar);
addParameter(p, 'vis', true, @islogical);
addParameter(p, 'figdir', fullfile('figures', mdlName), @ischar);
addParameter(p, 'mtfType', 'all', @(x) ismember(x, {'FP', 'LC', 'all'}));
addParameter(p, 'exclude_origin', true, @islogical);
addParameter(p, 'exclude_n_LC', 2, @isnumeric);  % Exclude models with more than this number of LCs
addParameter(p, 'LCGhostOnly', true, @islogical);  % Only include the slowest point on LC
addParameter(p, 'normalize', true, @islogical);
addParameter(p, 'demean', false, @islogical);
addParameter(p, 'Klist', 2:10, @isnumeric);
addParameter(p, 'K', 4, @isnumeric);
parse(p, varargin{:});

mtfFile = p.Results.mtfFile;
vis = p.Results.vis;
figdir = p.Results.figdir;
mtfType = p.Results.mtfType;
exclude_origin = p.Results.exclude_origin;
exclude_n_LC = p.Results.exclude_n_LC;
LCGhostOnly = p.Results.LCGhostOnly;
normalize = p.Results.normalize;
demean = p.Results.demean;
Klist = p.Results.Klist;
K = p.Results.K;

% infix for normalize
if normalize
    normstr = 'normalized';
else
    normstr = 'unscaled';
end

% infix for demeaning
if demean
    dmstr = 'demeaned';
else
    dmstr = 'uncentered';
end

if vis && ~exist(figdir, "dir")
    mkdir(figdir);
end


%% Load data and visualize distribution

load(mtfFile, 'allEq', 'allLCMotif');
try  % For compatibility with outdated data files
    size(allLCMotif);
catch
    load(mtfFile, 'allLC');
    allLCMotif = cellfun(@LC2Motif, allLC, 'UniformOutput', false);
    clear allLC
end

if vis && size(allEq, 2) > 1

    % Count number of motifs
    nFp = cellfun(@(x) size(x, 2), allEq);
    nLC = cellfun(@(x) size(x, 2), allLCMotif);
    nFp = categorical(string(nFp) + "FP");
    nLC = categorical(string(nLC) + "LC");
    nFpLC = nFp .* nLC;
    nFpLC = removecats(nFpLC);

    % Recode dynamics that are too rare as "others"
    pthres = 5e-3;
    [N, C] = histcounts(nFpLC);
    [N, idx] = sort(N, 'descend');
    C = C(idx);
    idx = find(N < pthres * numel(nFpLC), 1);
    nFpLC = mergecats(nFpLC, C(idx:end), 'others');
    C = [C(1:idx - 1) {'others'}];
    nFpLC = reordercats(nFpLC, C);

    [C, ~, ic] = unique(nFpLC);
    ic = reshape(ic, size(nFpLC));

    % Count
    nC = numel(C);
    N = zeros(nC);
    for i = 1:size(nFpLC, 1)
        N(ic(i, 1), ic(i, 2)) = N(ic(i, 1), ic(i, 2)) + 1;
    end
    N = N / size(nFpLC, 1);
    N1 = sum(N, 2);  % Total for session 1
    N2 = sum(N, 1);  % Total for session 2

    % Plot

    f = figure('Units', 'inches', 'Position', [0 0.5 9 6]);
    t = tiledlayout(nC + 1, nC + 1, "TileSpacing", "compact");

    % Distribution
    nexttile(1, [nC, nC]);
    myheatmap(N);
    colormap hot;
    cb = colorbar;
    cb.Layout.Tile = 'west';
    cb.FontSize = 12;
%     xticks([]);
    xticks(1:nC);
    xticklabels(C);
%     xlabel('Session 2');
    yticks(1:nC);
    yticklabels(C);
%     ytickangle(45);
    ylabel('Session 1');
    ax = gca();
    ax.YLabel.FontSize = 14;

    % Marginal for session 2
    cl = [min([N1(:); N2(:)]) max([N1(:); N2(:)])];
    nexttile(nC * (nC + 1) + 1, [1 nC]);
    myheatmap(N2);
    colormap(gca(), summer);
    clim(cl); colorbar off
%     xticks(1:nC);
%     xticklabels(C);
%     xlabel('Session 2')
    xticks([]);
    yticks(1);
    yticklabels("Total");

    % Marginal for session 1
    nexttile(nC + 1, [nC 1]);
    myheatmap(N1);
    colormap(gca(), summer);
    clim(cl);
    cb = colorbar;
    cb.Layout.Tile = 'east';
    cb.FontSize = 12;
    yticks([])
    xticks(1);
    xticklabels("Total");

    xlabel(t, 'Session 2');
%     ylabel(t, 'Session 1');
    t.XLabel.FontSize = 14;
%     t.YLabel.FontSize = 14;
    title(t, sprintf('Distribution of attractor type combinations across participants'), ...
        'FontSize', 16, 'FontWeight', 'bold');
%     subtitle(t, sprintf('(Attractor types with frequency < %g were grouped into "others")', pthres));
    PrintAsSeen(f, fullfile(figdir, ['motif_dist_' mdlName]), '-dsvg', '-vector');
    close(f);

end


%% Combine motifs

% Get names of motifs
nSubs = size(allEq, 1);
nSess = size(allEq, 2);
eqNames = cell(nSubs, nSess);
lcNames = cell(nSubs, nSess);
if nSess > 1
    for i = 1:nSubs
        for j = 1:nSess
            eqNames{i, j} = repmat(sprintf("S%dM%d-FP", i, j), 1, size(allEq{i, j}, 2));
            nLC = size(allLCMotif{i, j}, 2);
            lcNames{i, j} = cell(1, nLC);
            for k = 1:nLC
                lcNames{i, j}{k} = sprintf("S%dM%dLC%d-", i, j, k) + ...
                    ["vA" "vB" "vA" "vB"];  % A: major axis; B: minor axis
            end
        end
    end
else
    for i = 1:nSubs
        eqNames{i} = repmat(sprintf("S%d-FP", i), 1, size(allEq{i}, 2));
        nLC = size(allLCMotif{i}, 2);
        lcNames{i} = cell(1, nLC);
        for k = 1:nLC
            lcNames{i}{k} = sprintf("S%dLC%d-", i, k) + ...
                ["vA" "vB" "vA" "vB"];  % A: major axis; B: minor axis
        end
    end
end

% Exclusion
if exclude_origin
    for i = 1:numel(allEq)
        if isempty(allEq{i})
            continue;
        end
        idx = ~any(allEq{i});
        allEq{i}(:, idx) = [];
        eqNames{i}(:, idx) = [];
    end
end
if exclude_n_LC > 0
    for i = 1:numel(allLCMotif)
        if size(allLCMotif{i}, 2) > exclude_n_LC
            allLCMotif{i} = cell(1, 0);
            lcNames{i} = cell(1, 0);
            allEq{i} = [];
            eqNames{i} = [];
        end
    end
end
if LCGhostOnly
    for i = 1:numel(allLCMotif)
        if size(allLCMotif{i}, 2) == 1
            allLCMotif{i}{1} = allLCMotif{i}{1}(:, [1 3]);
            lcNames{i}{1} = lcNames{i}{1}([1 3]);
        elseif size(allLCMotif{i}, 2) > 1
            allLCMotif{i} = cellfun(@(x) x(:, 1), allLCMotif{i}, 'UniformOutput', false);
            lcNames{i} = cellfun(@(x) x(1), lcNames{i}, 'UniformOutput', false);
        end
    end
end

% Combine LCs for each model
for i = 1:numel(allLCMotif)
    allLCMotif{i} = [allLCMotif{i}{:}];
    lcNames{i} = [lcNames{i}{:}];
end

% Normalize the length of motifs
if normalize
    allEq = cellfun(@(x) x ./ vecnorm(x), allEq, 'UniformOutput', false);
    allLCMotif = cellfun(@(x) x ./ vecnorm(x), allLCMotif, 'UniformOutput', false);
end

% Demean
if demean
    allEq = cellfun(@(x) x - mean(x), allEq, 'UniformOutput', false);
    allLCMotif = cellfun(@(x) x - mean(x), allLCMotif, 'UniformOutput', false);
end

% Combine all motifs
allV = cellfun(@(x, y) [x y], allEq, allLCMotif, 'UniformOutput', false);
allVNames = cellfun(@(x, y) [x y], eqNames, lcNames, 'UniformOutput', false);


%% Select number of clusters

if vis && numel(Klist) > 1
    nk = numel(Klist);
    S = nan(nk, 3);
    semS = nan(nk, 3);

    for i = 1:nk
        [S(i, 1), semS(i, 1)] = ClusterInstabilityIndex([allV{:}]', Klist(i), [], ...
            60, 100, 'Distance', 'cosine', 'Replicates', 5);
        [S(i, 2), semS(i, 2)] = ClusterInstabilityIndex([allEq{:}]', Klist(i), [], ...
            60, 100, 'Distance', 'cosine', 'Replicates', 5);
        [S(i, 3), semS(i, 3)] = ClusterInstabilityIndex([allLCMotif{:}]', Klist(i), [], ...
            60, 100, 'Distance', 'cosine', 'Replicates', 5);
    end

    f = figure;
    errorbar(Klist, S, semS);
    xlabel('Number of Clusters');
    ylabel('Clustering Instability Index');
    legend({'Fixed points & ghosts', 'Only fixed points', 'Only ghosts'})
    title('Clustering Instability')
    PrintAsSeen(f, fullfile(figdir, [normstr '_' dmstr '_Klist_' mdlName]), '-dpng', '-r300');
end


%% Cluster motifs

switch mtfType
    case 'FP'
        V = allEq;
        VNames = eqNames;
    case 'LC'
        V = allLCMotif;
        VNames = lcNames;
    case 'all'
        V = allV;
        VNames = allVNames;
end

nParcels = size(V{find(cellfun(@(x) ~isempty(x), V), 1)}, 1);
atlas = GetYeoNetworks(nParcels);
parcelNames = atlas.Names(atlas.Net);
if vis
    [clIdx, C, f1, f2] = ClusterMotifs(V, VNames, K, 0, true, parcelNames, atlas.SortInd);
    title(f1.Children, sprintf("K-means clustering of %s motifs (%s)", mtfType, normstr))
    subtitle(f1.Children, mdlName, 'Interpreter', 'none');
    PrintAsSeen(f1, fullfile(figdir, ['K' num2str(K) '_' mtfType '_' normstr '_' dmstr '_motifs_' mdlName]), '-dpng', '-r300');
    PrintAsSeen(f2, fullfile(figdir, ['K' num2str(K) '_' mtfType '_' normstr '_' dmstr '_cluster_stats_' mdlName]), '-dpng', '-r300');
else
    [clIdx, C] = ClusterMotifs(V, VNames, K, 0, false);
end

% Reshape clIdx to (nSubs, nSess) cell of (1, nij) indices
nijmtf = cellfun(@(x) size(x, 2), V);
clIdx = mat2cell(clIdx', 1, nijmtf(:));
clIdx = reshape(clIdx, size(V));


%% Model parcel activation

if vis
    y = [V{:}];  % (nParcels, ?)
    parcelID = categorical(repmat(atlas.ParcelNames', 1, size(y, 2)));
    netID = categorical(repmat(atlas.Names(atlas.Net)', 1, size(y, 2)));
    clusterID = categorical(repmat([clIdx{:}], nParcels, 1));
    tbl = table(clusterID(:), netID(:), parcelID(:), y(:), 'VariableNames', ...
        {'cluster', 'net', 'parcel', 'y'});
    mdl = fitlme(tbl, "y ~ cluster * net + (cluster|net:parcel)");
    fixedMdl = anova(tbl, "y ~ cluster * net", "CategoricalFactors", {'cluster', 'net'}, ...
        'SumOfSquaresType', 'one');
    SS = stats(fixedMdl);
    SS = SS(:, "SumOfSquares");
    SS = [SS; {SS{"Error", 1} - mdl.SSE}];
    SS{"Error", 1} = mdl.SSE;
    SS.Row{end} = '(cluster|net:parcel)';
    SS = SS({'cluster', 'net', 'cluster:net', '(cluster|net:parcel)', 'Error', 'Total'}, :);
    SSv = [SS{1:end - 1, :}] / SS{"Total", 1} * 100;
%     randMdl = fitlme(tbl, "y ~ cluster + (cluster|net) + (cluster|net:parcel)");
%     reducedMdl = fitlme(tbl, "y ~ cluster + (cluster|net:parcel)");
%     compare(reducedMdl, mdl, "CheckNesting", true)

    figure('Units', 'inches', 'Position', [0 0.5 8 6]);
    bar(SSv);
    xticklabels(SS.Row(1:end - 1));
    ylabel('Explained sum of squares (%)');
    fontsize(gca, 14, 'points');
    title('Variation of attractors explained by each factor (type I SS)', 'FontSize', 18, ...
        'FontWeight', 'bold');
    subtitle("y ~ cluster + net + cluster:net + (cluster|net:parcel)", "FontSize", 16, "FontAngle", "italic");
    PrintAsSeen(fullfile(figdir, ['K' num2str(K) '_' mtfType '_' normstr '_' dmstr '_parcel_modeling_' mdlName]), '-dpdf', '-vector');
end


%% Visualize parcel activation similarity

if vis

    simMat = corrcoef([V{:}]');
    load(fullfile('data', 'atlas', 'Yeo_02_surf_dist.mat'), 'meanDist');
    atlas17 = GetYeoNetworks(nParcels, 17);
    atlas7 = GetYeoNetworks(nParcels, 7);
    Id17 = (atlas17.Net' == atlas17.Net);
    Id7 = (atlas7.Net' == atlas7.Net);

    % Some transformation
    Zr = atanh(simMat);
    S = -meanDist;

    % Regression
    tbl_full = table(S(:), Id7(:), Id17(:), Zr(:), 'VariableNames', {'S', 'Id7', 'Id17', 'Zr'});
    tbl_full = tbl_full(~isnan(tbl_full.S), :);
    tbl = tbl_full(:, {'S', 'Id17', 'Zr'});
    stepMdl = stepwiselm(tbl);
    aov = anova(stepMdl, 'component', '1');
    TSS = sum(aov.SumSq);
    aovId1 = anova(stepwiselm(tbl(:, [2 1 3])), 'component', '1');
    aovId1 = aovId1({'S', 'Id17', 'Id17:S', 'Error'}, :);
    aovSS3 = anova(stepMdl, 'component', '3');
    aovSS3 = aovSS3({'S', 'Id17', 'S:Id17', 'Error'}, :);

%     % Add inter-hemisphere distance to meanDist by epsilateral distance
%     % (only for visualization)
%     meanDist((end / 2 + 1):end, 1:(end / 2)) = meanDist(1:(end / 2), 1:(end / 2));
%     meanDist(1:(end / 2), (end / 2 + 1):end) = meanDist(1:(end / 2), 1:(end / 2));

    f = figure('Units', 'inches', 'Position', [0.5 0.5 12 9]);
    t = tiledlayout(2, 2, "TileSpacing", "loose");
    
    nexttile();
    MyNetMatYeo(Zr);
    colormap turbo
    cb = colorbar;
    cb.Label.String = 'atanh(\rho)';
    axis square
    title('Parcel activation similarity (Fisher Zr)', 'FontSize', 16, 'FontWeight', 'bold');

    nexttile();
    MyNetMatYeo(-meanDist);
    colormap turbo
    cb = colorbar;
    cb.Label.String = 'Distance (mm)';
    axis square
    title('Negative cortical distance (S = -D_{mesh})', 'FontSize', 16, 'FontWeight', 'bold');
%     subtitle(sprintf('SS(S | Intercept) = %.1f%%', aov.SumSq(1) / SS * 100), 'FontSize', 16, 'FontAngle', 'italic');

%     nexttile();
%     MyNetMatYeo(Id7);
%     colormap turbo
%     colorbar off
%     axis square
%     title('Same or different network (Id_7)', 'FontSize', 16, 'FontWeight', 'bold');
%     subtitle(sprintf('SS(Id_{7} | Intercept, S) = %.1f%%', aov.SumSq(2) / SS * 100), 'FontSize', 16, 'FontAngle', 'italic');

    nexttile();
    MyNetMatYeo(Id17);
    colormap turbo
%     colorbar off
    axis square
    title('Same or different network (Id_{17})', 'FontSize', 16, 'FontWeight', 'bold');
%     subtitle(sprintf('SS(Id_{17} | Intercept, S, Id_{7}) = %.1f%%', aov.SumSq(3) / SS * 100), 'FontSize', 16, 'FontAngle', 'italic');

    nexttile();
    SS = [aov.SumSq aovId1.SumSq aovSS3.SumSq];
    SS = SS / TSS * 100;
    bar(SS', 'stacked');
    hold on;
    for i = 1:size(SS, 2)
        tmp = [0; cumsum(SS(:, i))];
        for j = 1:size(SS, 1)
            if SS(j, i) < 1
                continue
            end
            text(i, mean(tmp(j:j+1)), sprintf("%.1f%%", SS(j, i)), ...
                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
        end
    end
    xticklabels({'Type I, S first', 'Type I, Id17 first', 'Type III'});
    ylabel('Percentage explained');
    legend({'S', 'Id17', 'Id17:S', 'Error'});
    fontsize(gca, 14, 'Points')
    title('Variation explained by each factor', 'FontSize', 16, 'FontWeight', 'bold')

%     nexttile(5, [1 2]);
%     histogram(simMat);
%     xlabel('Correlation (r)');
%     fontsize(gca, 14, 'points');
%     title('Distribution of activation similarity', 'FontSize', 16, 'FontWeight', 'bold')
% 
%     nexttile(7, [1 2]);
%     histogram(S);
%     xlabel('Negative cortical surface distance (mm)');
%     fontsize(gca, 14, 'points');
%     title(['Distribution of negative cortical distance' tmpStr], 'FontSize', 16, 'FontWeight', 'bold')

    title(t, "Similarity of parcel activation across attractors", 'FontSize', 18, 'FontWeight', 'bold');
    subtitle(t, "Zr ~ 1 + S + Id_{17}, stepwise", 'FontSize', 16, 'FontAngle', 'italic');

    PrintAsSeen(f, fullfile(figdir, ['K' num2str(K) '_' mtfType '_' normstr '_' dmstr '_parcel_similarity_' mdlName]), '-dpdf', '-vector');
end


%% Visualize network activation on same axis

if vis
    f = figure('Units', 'inches', 'Position', [0 0.5 20 10.0208]);
    X = [allV{:}];  % (nParcels, N)
    cIdx = [clIdx{:}];  % (1, N)
    nNet = max(atlas.Net);
    nIdx = atlas.Net';
    c = distinguishable_colors(nNet);
    for i = 1:K
        nexttile();
        hold on
        for j = 1:nNet
            tmpDat = X(nIdx == j, cIdx == i);
            tb = boxchart(ones(numel(tmpDat), 1) * j, tmpDat(:), 'Notch', 'on');
            tb.BoxFaceColor = c(j, :);
            tb.MarkerColor = c(j, :);
            tb.WhiskerLineColor = c(j, :);
        end
        legend off;
        xticks(1:nNet);
        xticklabels(atlas.Names);
        yline(0, '--k');
        fontsize(gca, 14, "points");
        title("Cluster " + i, 'FontSize', 16, 'FontWeight', 'bold');
    end
    lgd = legend(atlas.Names);
    lgd.Layout.Tile = 'east';
    lgd.FontSize = 14;
    title(gcf().Children, sprintf('Network activation in each cluster'), ...
        'FontSize', 18, 'FontWeight', 'bold');
    subtitle(gcf().Children, mdlName, 'Interpreter', 'none', ...
        'FontSize', 16, 'FontAngle', 'italic');
    PrintAsSeen(fullfile(figdir, ...
        ['K' num2str(K) '_' mtfType '_' normstr '_' dmstr '_networks_origOrder_' mdlName]), '-dpng', '-r300');
end


%% Visualizing network activation

if vis
    X = [allV{:}];  % (nParcels, N)
    Y = [clIdx{:}];  % (1, N)
    nid = atlas.Net';
    nNet = max(nid);
    netAct = nan(nNet, K);
    for i = 1:nNet
        for j = 1:K
            netAct(i, j) = mean(X(nid == i, Y == j), 'all');
        end
    end

    c = distinguishable_colors(nNet);
    f = figure('Units', 'inches', 'Position', [0 0.5 20 10.0208]);
    for i = 1:K
        nexttile();
        hold on
        [~, newNetOrder] = sort(netAct(:, i), 'descend');
        for j = 1:nNet
            tmpDat = X(nid == newNetOrder(j), Y == i);
            tb = boxchart(ones(numel(tmpDat), 1) * j, tmpDat(:), 'Notch', 'on');
            tb.BoxFaceColor = c(newNetOrder(j), :);
            tb.MarkerColor = c(newNetOrder(j), :);
            tb.WhiskerLineColor = c(newNetOrder(j), :);
        end
        legend off;
        xticks(1:nNet);
        xticklabels(atlas.Names(newNetOrder));
        yline(0, '--k');
        fontsize(gca, 14, "points");
        title("Cluster " + i, 'FontSize', 16, 'FontWeight', 'bold');
    end
    lgd = legend(atlas.Names(newNetOrder));
    lgd.Layout.Tile = 'east';
    lgd.FontSize = 14;
    title(gcf().Children, sprintf('Network activation in each cluster'), ...
        'FontSize', 18, 'FontWeight', 'bold');
    subtitle(gcf().Children, mdlName, 'Interpreter', 'none', ...
        'FontSize', 16, 'FontAngle', 'italic');
    PrintAsSeen(fullfile(figdir, ...
        ['K' num2str(K) '_' mtfType '_' normstr '_' dmstr '_networks_' mdlName]), '-dpng', '-r300');

    % Same plot but only look at cluster centroids
    c = distinguishable_colors(nNet);
    f = figure('Units', 'inches', 'Position', [0 0.5 20 10.0208]);
    for i = 1:K
        nexttile();
        hold on
        [~, newNetOrder] = sort(netAct(:, i), 'descend');
        for j = 1:nNet
            tmpDat = C(i, nid == newNetOrder(j))';
            tb = boxchart(ones(numel(tmpDat), 1) * j, tmpDat(:), 'Notch', 'on');
            tb.BoxFaceColor = c(newNetOrder(j), :);
            tb.MarkerColor = c(newNetOrder(j), :);
            tb.WhiskerLineColor = c(newNetOrder(j), :);
        end
        legend off;
        xticks(1:nNet);
        xticklabels(atlas.Names(newNetOrder));
        yline(0, '--k');
        fontsize(gca, 14, "points");
        title("Centroid " + i, 'FontSize', 16, 'FontWeight', 'bold');
    end
    lgd = legend(atlas.Names(newNetOrder));
    lgd.Layout.Tile = 'east';
    lgd.FontSize = 14;
    title(gcf().Children, sprintf('Network activation in each cluster centroid (across parcels)'), ...
        'FontSize', 18, 'FontWeight', 'bold');
    subtitle(gcf().Children, mdlName, 'Interpreter', 'none', ...
        'FontSize', 16, 'FontAngle', 'italic');
    PrintAsSeen(fullfile(figdir, ...
        ['K' num2str(K) '_' mtfType '_' normstr '_' dmstr '_centroid_networks_' mdlName]), '-dpng', '-r300');

end


%% Close the figures;

close all

end


%% Utilities

function ax = myheatmap(N)
    hold on
    imagesc(N);
    set(gca, 'YDir','reverse');
    cl = clim();
    for i = 1:size(N, 1)
        for j = 1:size(N, 2)
            if N(i, j) > 0.2 * cl(1) + 0.8 * cl(2)
                c = 'k';
            else
                c = 'w';
            end
            text(j, i, sprintf('%.4g', N(i, j)), "FontSize", 12, "Color", c, ...
                "HorizontalAlignment", "center", "VerticalAlignment", "middle");
        end
    end
    axis tight;
    ax = gca; % get the current axes
    ax.XAxis.FontSize = 12; % change the font size of the X tick labels
    ax.YAxis.FontSize = 12; % change the font size of the Y tick labels
end