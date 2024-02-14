function [motifs, centroids, dist, Ct, Cp, subjType] = EqLC2Motifs(allEq, allLC, K, Distance)
%EQLC2MOTIFS Categorize subject motifs
%
%   Inputs:
%     - |allEq|, |allLC|: from attractor analysis.
%     - |K|: number of different motifs (w\o reflections), default 2
%     - |Distance|: metric for |kmeans()|. Default is cosine.
%
%   Output |motifs| is a (nParcels, K, nSubs, nSess) array. For each model,
%   each motif, the attractor/LCPC of this model that is most similar to
%   the centroid is selected REGARDLESS OF WHETHER IT WAS CLUSTERED INTO
%   THIS GROUP. Therefore, for people that do not have a motif in each
%   cluster, the selected vectors for a cluster may not belong to this
%   cluster.
%
%   |centroids| is a (K, nParcels) matrix, the centroids of each cluster.
%
%   |dist| is a (K, nSubs, nSess) array, the distance of |motifs| to
%   the associated centroid.
%
%   |Ct| is a (K, nParcels) matrix, the t-statistics for the activation of
%   each parcel within each cluster (ONLY based on the vectors that are
%   really in the cluster). |Cp| is the associated uncorrected p-value.

if nargin < 4 || isempty(Distance)
    Distance = 'cosine';
end
if nargin < 3 || isempty(K)
    K = 2;
end
nSubs = size(allEq, 1);
nSess = size(allEq, 2);

% Limit cycle PC1&2 scaled to maximum projection length, w/ reflections
LCPCs = cell(nSubs, nSess);
for i = 1:numel(allLC)
    if ~isempty(allLC{i})
        curr = [allLC{i}{:}];
        [coeff, scores] = pca(curr', 'NumComponents', 2);
        LCPC = coeff .* max(abs(scores));
        LCPCs{i} = [LCPC -LCPC];  % Include reflections
    end
end

% All vectors, including reflections
allV = cell(nSubs, nSess);
for i = 1:numel(allEq)
    allV{i} = [allEq{i} LCPCs{i}];  % allEq already includes reflections
end

% Clustering
tmp = [allV{:}]';
[clIdx, C, ~, D] = kmeans(tmp, K * 2, 'Replicates', 10, 'Distance', Distance);

% Calculate a t-value to show consistency of each parcel
Ct = nan(size(C));
Cp = nan(size(C));
for k = 1:K * 2
    [~, Cp(k, :), ~, stat] = ttest(tmp(clIdx == k, :));
    Ct(k, :) = stat.tstat;
end

% Identify motifs with mostly negative sign
corrMat = corrcoef(C');
[~, pairID] = min(corrMat);  % (i, pairID(i)) is the same motif with different sign
pairs = unique([min(1:K * 2, pairID); max(1:K * 2, pairID)]', "rows");
motifCID = nan(K, 1);
for i = 1:K
    if mean(sign(C(pairs(i, 1), :))) <= 0
        motifCID(i) = pairs(i, 1);
    else
        motifCID(i) = pairs(i, 2);
    end
end

% Relabel according to popularity
clSizes = sum(clIdx == motifCID');
[~, tmpIdx] = sort(clSizes, 'descend');  % e.g., tmpIdx = [2 4 3 1] (cluster2 is most popular)
motifCID = motifCID(tmpIdx);

% Recode the clusters
centroids = C(motifCID, :);
trueD = min(D(:, motifCID), D(:, pairID(motifCID)));  % Distance, ignoring reflection
[~, trueIdx] = min(trueD, [], 2);   % "True" cluster label, ignoring reflection
D = D(:, motifCID);
Ct = Ct(motifCID, :);
Cp = Cp(motifCID, :);

% Transform to cell due to diferent number of vectors in each cell
nV = cellfun(@(x) size(x, 2), allV);
D = mat2cell(D, nV(:), K);
D = reshape(D, nSubs, nSess);

% Pick the most similar vector
nParcels = size(centroids, 2);
motifs = nan(nParcels, K, nSubs, nSess);
dist = nan(K, nSubs, nSess);
for i = 1:nSubs
    for j = 1:nSess
        [dist(:, i, j), tmp] = min(D{i, j});  % D{i, j} is of size (nV(i, j), K)
        motifs(:, :, i, j) = allV{i, j}(:, tmp);
    end
end

% Calculate subject type
trueIdx = mat2cell(trueIdx, nV(:), 1);
trueIdx = reshape(trueIdx, nSubs, nSess);
subjType = cellfun(@(x) num2str(unique(x'), 'M%d'), trueIdx, 'UniformOutput', false);
subjType = categorical(subjType);

end