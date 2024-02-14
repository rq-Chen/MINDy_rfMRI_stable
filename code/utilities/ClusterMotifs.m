function [clIdx, centroids, f1, f2] = ClusterMotifs(allV, allVNames, nCl, sgn, vis, parcelNames, sortInd)
%CLUSTERMOTIFS Cluster the motifs (fixed points, limit cycle features)
%
% Note: since the vectors have unit length, it makes sense to define their
% distance by their geodesic distance (on a 100-dimensional hyperball),
% which is equivalent to their angular distance, or their cosine distance.
%
% Inputs:
%   - allV: (nParcels, N), a set of N unit vectors, or cell array of such
%   sets of vectors
%   - allVNames: a string array of the same length as the total number of
%   vectors, or a cell array of string arrays that has a total number of
%   this many strings. Make sure to use strings instead of char arrays!
%   - nCl: number of clusters (default 3)
%   - sgn: force the majority of elements to have this sign (0 for no flipping)
%   - vis: do visualization or not (dafault true)
%   - parcelNames: names of the parcels (pre-sorting)
%   - sortInd: reorder the parcels when do plot
%
% Outputs:
%   - clIdx: (nMotifs, 1) cluster indices
%   - centroids: (nCl, nParcels) cluster centroids
%   - f1: handle of the first figure (coefficients and centroids)
%   - f2: handle of the second figure (distribution of clusters across
%   motif types and sessions)

% Data and labels
if iscell(allV)
    allV = [allV{:}];
end
if nargin < 2 || isempty(allVNames)
    allVNames = string(1:size(allV, 2));
elseif iscell(allVNames)
    allVNames = [allVNames{:}];
end
assert(length(allVNames) == size(allV, 2), ...
    "Incompatible number of vectors and labels. Make sure the labels are strings, not char vectors!")
nParcels = size(allV, 1);
if iscolumn(allVNames); allVNames = allVNames'; end

% Hyperparameters
if nargin < 4 || isempty(sgn)
    sgn = -1;  % 1 for mostly positive, -1 for mostly negative, 0 for no flipping
end
if nargin < 3 || isempty(nCl)
    nCl = 3;  % Number of clusters
end

% Visualization
if nargin < 5 || isempty(vis)
    vis = true;
end
if vis
    if nargin < 6 || isempty(parcelNames)
        parcelNames = string(1:nParcels);
    end
    if nargin < 7 || isempty(sortInd)
        sortInd = 1:nParcels;
    end
end

% Flip sign
if sgn
    tmp = sign(mean(sign(allV)));
    tmp(tmp == 0) = 1;
    allV = allV .* (sgn * tmp);
end

% Clustering
[clIdx, centroids] = kmeans(allV', nCl, 'Replicates', 10, 'Distance', 'cosine');

% Relabel the clusters by popularity
clSizes = sum(clIdx == (1:nCl));
[clSizes, tmpIdx] = sort(clSizes, 'descend');  % e.g., tmpIdx = [2 4 3 1] (cluster2 is most popular)
centroids = centroids(tmpIdx, :);
[~, tmpIdx] = sort(tmpIdx);  % e.g., tmpIdx = [4 1 3 2]
clIdx = tmpIdx(clIdx);  % cluster2 -> 1, cluster 4 -> 2, cluster 3 -> 3, etc.
clIdx = clIdx';  % (nMotifs, 1) like the output of kmeans

% Sort the data with the clusters
[~, idx] = sort(clIdx);

% Visualization
if ~vis
    return
end
basefont = 14;
nmtfcols = 2;
nsurfcols = 2 - (nCl <= 2);
ncols = nmtfcols + nsurfcols;
nrows = ceil(nCl / nsurfcols);

f1 = figure('Units', 'inches', 'Position', [0 0.5 20 10.0208]);
t = tiledlayout(nrows, ncols, 'TileSpacing', 'tight', 'Padding', 'tight');
colormap parula
nexttile(1, [nrows nmtfcols])
hold on;
imagesc(allV(sortInd, idx));
% clim([-1 1] * max(abs(clim())));  % Symmetric color scale
colorbar;

% Cluster and network boundaries
[unq, ~, ic] = unique(parcelNames(sortInd), 'stable');
netIStart = arrayfun(@(x) find(ic == x, 1, 'first'), 1:numel(unq));
yline(netIStart(2:end) - 0.5, 'k', 'LineWidth', 1);
xline(cumsum(clSizes(1:end - 1)) + 0.5, 'k', 'LineWidth', 2);
axis tight

% Ticks
if numel(allVNames) < 100
    xticks(1:size(allV, 2));
    xticklabels(allVNames(idx));
    xtickangle(90);
else  % Add a label for each cluster at the mean index
    iStart = [0, cumsum(clSizes)] + 1;
    xticks((iStart(1:end - 1) + iStart(2:end) - 1) / 2);
    xticklabels("Cluster " + string(1:nCl));
end
if numel(parcelNames) < 100
    yticks(1:nParcels);
    yticklabels(parcelNames(sortInd));
else  % Add a label for each network at the mean index
    iStart = [netIStart, nParcels + 1];
    yticks((iStart(1:end - 1) + iStart(2:end) - 1) / 2);
    yticklabels(unq);
end
fontsize(gca, basefont, "points");
title('Coefficients', 'FontSize', basefont + 2);

% Surface plots
% cl = [-1 1] * max(abs(centroids), [], 'all');
cl = [min(centroids, [], "all") max(centroids, [], "all")];
for i = 1:nCl
    nexttile(tilenum(t, ceil(i / nsurfcols), mod(i - 1, nsurfcols) + 1 + nmtfcols));
    PlotYeoSurface(centroids(i, :));
    clim(cl);
    fontsize(gca, basefont, "points");
    title(['Cluster ' num2str(i)], 'FontSize', basefont + 2)
end
cb = colorbar; cb.Layout.Tile = 'east';
title(t, 'K-means clustering of the motifs', 'FontSize', basefont + 4, 'FontWeight', 'bold')
subtitle(t, sprintf('%d parcel models', nParcels), 'FontSize', basefont + 2, 'FontAngle', 'italic');

% figure;
% score = myPCA(allV', 2);
% % score = mdscale(squareform(pdist(allV)), 2);
% gscatter(score(:, 1), score(:, 2), clIdx);
% xlabel('PC1'); ylabel('PC2');
% title('Dynamical motifs in PC space')

% The distribution of different vectors
words = split(allVNames', "-");
words = unique(words(:, end))';
countMat = nan(nCl, length(words));
vType = nan(numel(allVNames), numel(words));
for i = 1:numel(words)
    vType(:, i) = contains(allVNames, words(i));
end
for i = 1:nCl
    countMat(i, :) = sum(vType(clIdx == i, :));
end
f2 = figure;
heatmap(countMat, "XDisplayLabels", words, ...
    "YDisplayLabels", "Cluster " + string(1:nCl), ...
    "Title", "Distribution of motifs within each cluster")
end