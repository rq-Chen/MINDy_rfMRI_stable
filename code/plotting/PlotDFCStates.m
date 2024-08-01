function [f] = PlotDFCStates(C1, varargin)
%PLOTDFCSTATES Plot the DFC states and transition probabilities
%
%   PlotDFCStates(C1);  % For one dataset
%   PlotDFCStates(C1, C2);  % For two datasets
%
% Inputs are outputs of `dFC_state_analysis` function. The number
% of clusters must be 1 to maxK.

% Check the number of input arguments
if nargin == 1
    isCmp = false;
elseif nargin == 2
    isCmp = true;
    C2 = varargin{1};
else
    error('The number of input arguments must be 1 or 2.');
end

% Reorder the states of C2 to align to C1
if isCmp
    for i = 2:numel(C2)
        sortIdx = reorderStates(C1{i}, C2{i});
        C2{i} = C2{i}(:, :, sortIdx);
    end
end

% Combine C1 and C2
if isCmp
    for i = 1:numel(C1)
        for j = 1:size(C1{i}, 3)
            C1{i}(:, :, j) = tril(C1{i}(:, :, j), -1) + triu(C2{i}(:, :, j), 1);
        end
    end
end

maxK = numel(C1);
nRows = maxK;
nCols = maxK;
f = figure('Units', 'inches', 'Position', [0 0 12 15]);
tlo = tiledlayout(nRows, nCols, "TileSpacing", "tight", "Padding", "compact");
cl = [-0.4, 0.8];

% Plot cluster centers
for i = 1:maxK
    for j = 1:i
        nexttile(tilenum(tlo, i, j));
        ShowXText = 'none';
        ShowYText = 'none';
        W = C1{i}(:, :, j);
        MyNetMatYeo(W + diag(nan(size(W, 1), 1)), [], ShowXText, ShowYText, false);
        clim(cl);
        colorbar off;
        axis square;
        if ~isCmp
            title(sprintf('State %d of %d', j, i));
        else
            title(sprintf('State %d of %d, \\rho = %.2f', j, i, ...
                FC_similarity(W, W')));
        end
    end
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Visible = 'on';
% title(tlo, 'DFC states')

end


%% Helper functions

function sortIdx = reorderStates(C1, C2, varargin)
%REORDERSTATES Reorder the states of C2 to align to C1
%
%   sortIdx = REORDERSTATES(C1, C2, varargin)
%
% Inputs:
%   C1: (nParcels, nParcels, K1) matrix
%   C2: (nParcels, nParcels, K2) matrix
%   varargin: Additional arguments for pdist2 function
%
% Outputs:
%   sortIdx: (K2, 1) vector. Reorder C2 by C2(:, :, sortIdx)
%     and TPM2 by TPM2(sortIdx, sortIdx).

% Convert matrices to vectors
C1 = Matrix2Vec(C1)';  % From dFCwalk toolbox
C2 = Matrix2Vec(C2)';
distMat = pdist2(C2, C1, varargin{:});  % (K2, K1) matrix
[tmpIdx, ~] = assignmentoptimal(distMat);  % Hungarian algorithm
tmpIdx(tmpIdx == 0) = inf;
[~, sortIdx] = sort(tmpIdx);

end


function sim = FC_similarity(FC1, FC2)
    sim = corr(Matrix2Vec(FC1), Matrix2Vec(FC2));
end