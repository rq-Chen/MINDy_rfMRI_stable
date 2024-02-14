function [] = PlotLimitSet(lmst, mdl, dim, maxUSDim)
%PLOTLIMITSET Plot the limitset of one model
%   Inputs:
%   - |lmst|: should be the content of one cell in |allLmSt| returned
%     by |CategorizeLimitSets|
%   - |mdl|: the corresponding model
%   - |dim|: visualization dimensionality (default 3)
%
%   Currently, for limit cycles the script will plot the limit cycles, but
%   for other models it will only plot the fixed points and simulation
%   starting around the origin.
if nargin < 4
    maxUSDim = 1;
end
if nargin < 3
    dim = 3;
end
if nargin < 2
    mdl = [];
end

% Perturbation parameters
nSims = 120;
ptbStd = .05;
T = 2000;

currHold = ishold; hold on;

if iscell(lmst)
    nLC = length(lmst);
    c = distinguishable_colors(nLC);
    allT = cellfun(@(x) size(x, 2), lmst);
    cumT = [0 cumsum(allT)];
    [~, score] = pca([lmst{:}]', 'NumComponents', dim);
    for i = 1:nLC
        if dim == 2
            plot(score(cumT(i) + 1:cumT(i + 1), 1), score(cumT(i) + 1:cumT(i + 1), 2), ...
                'Color', c(i, :));
        elseif dim == 3
            plot3(score(cumT(i) + 1:cumT(i + 1), 1), score(cumT(i) + 1:cumT(i + 1), 2), ...
                score(cumT(i) + 1:cumT(i + 1), 3), 'Color', c(i, :));            
        end
    end
elseif isstruct(lmst)
    adjMat = cellfun(@(x) ~isempty(x), lmst.E);
    % Include the origin but remove the non-seperatices
    fpIdx = (lmst.VUSDim == 0 | ...
        (lmst.VUSDim <= maxUSDim & any(adjMat, 2)) | ...
        vecnorm(lmst.V, 2, 2) < 1e-6);
    nFp = sum(fpIdx);
    nParcels = size(lmst.V, 2);

    c = distinguishable_colors(nFp);
    sims = MINDyInt_00(mdl, randn(nParcels, nSims) * ptbStd, 1, 1, 0, T);
    [score, coeff, ~, mu] = myPCA(sims, dim);
    FpScore = (lmst.V(fpIdx, :) - mu) * coeff;
    tmpDist = pdist2(squeeze(score(end, :, :)), FpScore);
    [~, tmpIdx] = min(tmpDist, [], 2);
    for i = 1:nSims
        if dim == 2
            plot(score(:, i, 1), score(:, i, 2), 'Color', c(tmpIdx(i), :));
        elseif dim == 3
            plot3(score(:, i, 1), score(:, i, 2), score(:, i, 3), 'Color', c(tmpIdx(i), :));
        end
    end
    FpDim = [lmst.VUSDim(fpIdx)];
    scatter(FpScore(:, 1), FpScore(:, 2), 200, c, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 2);
    text(FpScore(:, 1), FpScore(:, 2), num2str(FpDim), ...
        "FontSize", 32, 'BackgroundColor', [1 1 1 .9]);
end

if dim == 3; view([-37.5, 30]); end
if ~currHold; hold off; end

end

