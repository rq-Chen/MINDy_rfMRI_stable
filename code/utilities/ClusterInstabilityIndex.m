function [S, semS, SA, SR] = ClusterInstabilityIndex(X, k, f, nDiv, nRnd, varargin)
%CLUSTERINSTABILITYINDEX Instability index for clustering algorithms
%   Clustering instability index in (Lange et al., 2004) as used in
%   (Yeo et al., 2011) network parcellation.
%
%   Inputs:
%     - X: data, (nObservation, nFeatures) matrix
%     - k: number of clusters
%     - f: a function handle: (idx, fc) = f(X, k, varargin{:})
%       - idx: cluster index, (nObservation, 1) vector
%       - fc: a function handle: idx_new = fc(X_new) (the classifier)
%       By default, |f| is the kmeans algorithm.
%     - nDiv: number of random partitions (default is 30)
%     - nRnd: number of (pairs of) random labeling to compute S(R_k)
%       (default is 30 pairs).
%
%   Outputs:
%     - S: instability index, normalized by that of random labels
%     - semS: standard error of mean for S (over all partitions, 
%       assuming that SR is fixed at its mean)
%     - SA: instability index for the given f, (nDiv, 1) vector
%     - SR: instability index for random labeling, (nRnd, 1) vector
%
%   This script uses the Hungarian algorithm implemented in Markus
%   Buehren's |assignmentoptimal.m| function.

if nargin < 5 || isempty(nRnd)
    nRnd = 30;
end
if nargin < 4 || isempty(nDiv)
    nDiv = 30;
end
if nargin < 3 || isempty(f)
    f = @f_kmeans;
end
N = size(X, 1);

SA = nan(nDiv, 1);
for i = 1:nDiv
    permIdx = randperm(N, floor(N / 2));
    cmpIdx = setdiff(1:N, permIdx);
    [~, fc] = f(X(permIdx, :), k, varargin{:});
    idx = f(X(cmpIdx, :), k, varargin{:});
    pred = fc(X(cmpIdx, :));
    SA(i) = MisclassificationCost(idx, pred, k);
end

SR = nan(nRnd, 1);
for i = 1:nRnd
    idx = randi(k, size(cmpIdx));
    pred = randi(k, size(cmpIdx));
    SR(i) = MisclassificationCost(idx, pred, k);
end

S = mean(SA) / mean(SR);
semS = std(SA) / mean(SR) / sqrt(nDiv);

end


%% Misclassification cost

function cost = MisclassificationCost(idx, pred, k)

contingency = zeros(k);
n = numel(idx);
for j = 1:n
    contingency(idx(j), pred(j)) = contingency(idx(j), pred(j)) + 1;
end
costMat = sum(contingency) - contingency;
[~, cost] = assignmentoptimal(costMat);
cost = cost / n;

end


%% K-means function

function [idx, fc] = f_kmeans(X, k, varargin)
%F_KMEANS The utility function for kmeans algorithm

[idx, C] = kmeans(X, k, varargin{:});
fc = @localfunc;

tmp = find(cellfun(@(x) (ischar(x) || isstring(x)) && strcmpi(x, 'Distance'), varargin));
if tmp
    distance_type = varargin{tmp + 1};
else
    distance_type = 'euclidean';
end

    function idxnew = localfunc(xnew)
        dist = pdist2(C, xnew, distance_type);
        [~, idxnew] = min(dist);
        idxnew = idxnew';
    end

end



