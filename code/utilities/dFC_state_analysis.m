function [C, TPM, clIdx] = dFC_state_analysis(dat, allK, W, lag)
%DFC_STATE_ANALYSIS calculates dynamic functional connectivity states
%
%   [C, TPM, clIdx] = dFC_state_analysis(dat, allK, W, lag)
%
%   Inputs:
%       dat: (nParcels, nTRs) matrix or cell array of matrices
%       K: (nK, 1) SORTED array, number of states, default is 2:8
%       W: window size in seconds, default is 60s
%       lag: lag in seconds, default is 10s
%
%   Outputs:
%       C: (nK, 1) cell array of (nParcels, nParcels, K) matrices, cluster centroids
%       TPM: (nK, 1) cell array of (K, K) matrices, transition probability matrices
%       clIdx: (nK, 1) cell array of [(1, nTRs) cluster indices or cell array of them]
%
%   The FC matrices were Fisher z-transformed before clustering and the outputs were
%   also Fisher z-transformed.

if ~iscell(dat)
    dat = {dat};
end
if nargin < 2 || isempty(allK)
    allK = 2:8;
end
nK = numel(allK);
if nargin < 3 || isempty(W)
    W = 60;
end
if nargin < 4 || isempty(lag)
    lag = 10;
end

% Convert W and lag to TRs
TR = 0.72;
W = round(W / TR);
lag = round(lag / TR);

% Calculate dFC
dFC = cell(size(dat));
nW = nan(size(dat));
for i = 1:numel(dFC)
    dFC{i} = TS2dFCstream(dat{i}', W, lag, '2D');
    nW(i) = size(dFC{i}, 2);
end

% Fisher z-transform
dFC = cellfun(@(x) atanh(x), dFC, 'UniformOutput', false);

% Clustering
C = cell(nK, 1);
TPM = cell(nK, 1);
clIdx = cell(nK, 1);
for i = 1:nK
    if allK(i) == 1
        C{i} = mean([dFC{:}], 2)';
        TPM{i} = eye(1);
        continue
    end
    [idx, tmpC] = kmeans([dFC{:}]', allK(i), 'Replicates', 3);
    tmpTPM = idx2TPM(idx, allK(i), nW);

    % Align the centroids for allK(i - 1) and allK(i)
    if i > 1
        distMat = pdist2(tmpC, C{i - 1});  % (allK(i), allK(i - 1))
        [tmpIdx, ~] = assignmentoptimal(distMat);  % Hungarian algorithm
        origIdx = 1:allK(i);
        tmpIdx(tmpIdx == 0) = origIdx(tmpIdx == 0);
        [~, sortIdx] = sort(tmpIdx);
        tmpC = tmpC(sortIdx, :);
        tmpTPM = tmpTPM(sortIdx, sortIdx);
    end

    C{i} = tmpC;
    TPM{i} = tmpTPM;
    newIdx = tmpIdx(idx)';
    newIdx = reshape(mat2cell(newIdx, 1, nW(:)'), size(nW));
    clIdx{i} = newIdx;
end

% Convert to matrices
for i = 1:nK
    C{i} = mat2cell(C{i}', size(C{i}, 2), ones(1, size(C{i}, 1)));
    C{i} = cellfun(@Vec2Mat, C{i}, 'UniformOutput', false);
    C{i} = cat(3, C{i}{:});
end


end


%% Helper functions

function TPM = idx2TPM(idx, K, nW)
%IDX2TPM converts cluster indices to transition probability matrix
    idx = mat2cell(idx, nW(:));
    TPM = zeros(K, K);
    for i = 1:numel(idx)
        tmp = idx{i};
        for j = 1:numel(tmp)-1
            TPM(tmp(j), tmp(j+1)) = TPM(tmp(j), tmp(j+1)) + 1;
        end
    end
    TPM = TPM ./ sum(TPM, 2);
end


function A = Vec2Mat(V)
    A = squareform(V);
    A = A - diag(nan(size(A, 1), 1));
end