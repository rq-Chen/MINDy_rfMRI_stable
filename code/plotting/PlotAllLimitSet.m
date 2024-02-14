%% Show the graphs of limit sets
function clIdx = PlotAllLimitSet(nICs, allLmSt, mdlName, maxUSDim)

if nargin < 4
    maxUSDim = 1;
end
if nargin < 3
    mdlName = 'FinVal100';
end

% Sizes
nICs = permute(nICs, [1 3 2]);
[maxFpOnIC, nSess, nSub] = size(nICs);

% Overview of all model's invariant cycles
nICs = reshape(nICs, maxFpOnIC, []);
mdlLabel = "S" + string(1:nSub) + ("M" + string((1:nSess)'));
mdlLabel = mdlLabel(:);
cl = clustergram(nICs, 'ColumnLabels', mdlLabel, 'RowLabels', string([0 2:maxFpOnIC]) + "FPs", ...
    'Annotate', true, 'AnnotPrecision', 1, 'DisplayRange', 6, ...
    'Cluster', 'row', 'ColumnPDist', 'hamming');
addTitle(cl, sprintf('Subclasses of models based on number & types of invariant cycles, Model = %s', mdlName));

% Show representative plots for each cluster

clIdx = clusterdata(nICs', 'distance', 'hamming', 'criterion', 'distance', 'cutoff', 0.01);
nCl = max(clIdx);

% Re-label clusters from most to least popular
nMdlInCl = sum(clIdx == 1:nCl)';
[nMdlInCl, tmpIdx] = sort(nMdlInCl, 'descend');
[~, tmptmpIdx] = sort(tmpIdx);  % tmpIdx(tmptmpIdx(i)) = i;
clIdx = tmptmpIdx(clIdx);

ny = 2; nx = ceil(nCl / ny);
f = figure();
f.WindowState = 'maximized';
t = tiledlayout(ny, nx);
for icl = 1:nCl
    nexttile;
    tmp = find(clIdx == icl, 1);
    curr = allLmSt{ceil(tmp / nSess), 1 + mod(tmp - 1, nSess)};
    if iscell(curr)  % Limit Cycles
        nCy = numel(curr);
        if mod(nCy, 2)
            viscircles([0 0], 1, 'Color', 'k');
        end
        if nCy > 1
            ang = linspace(0, 2 * pi, nCy - mod(nCy, 2) + 1)';
            ang = ang(1:end - 1);
            viscircles(2 * [cos(ang) sin(ang)], 0.5, 'Color', 'k');
        end
        axis equal; xticks([]); yticks([]);
    else  % Fixed points

        % The nodes to plot
        adjMat = cellfun(@(x)~isempty(x), curr.E);
        idx = (any(adjMat | adjMat')' | curr.VUSDim == 0 | vecnorm(curr.V, 2, 2) < 1e-6); 

        c = distinguishable_colors(maxUSDim + 2);
        tmp = c(1, :); c(1, :) = c(2, :); c(2, :) = tmp;  % Use red for attractors
        cIdx = curr.VUSDim + 1;
        cIdx(cIdx > maxUSDim + 2) = maxUSDim + 2;
        c = c(cIdx, :);
        l = curr.VUSDim();
        
        adjMat = adjMat(idx, idx);
        c = c(idx, :);
        l = l(idx);

        p = plot(digraph(adjMat), 'NodeColor', c, ...
            'NodeLabel', l, 'Layout', 'force');
        [~, origIdx] = min(vecnorm(curr.V(idx, :), 2, 2));
        p.XData(origIdx) = mean(p.XData([1:origIdx - 1 origIdx + 1:length(p.XData)]));
        p.YData(origIdx) = mean(p.YData([1:origIdx - 1 origIdx + 1:length(p.YData)]));
        
%         % Compute a coordinate for the nodes
%         [~, coeff] = myPCA([curr.E{:}], 3);
%         if isempty(coeff)
%             coeff = pca(curr.V(idx, :), 'NumComponents', 3);
%             if size(coeff, 2) < 3
%                 coeff = [coeff zeros(size(coeff, 1), 1)];
%             end
%         end
%         scores = curr.V(idx, :) * coeff;
% 
%         plot(digraph(adjMat), 'NodeColor', c, ...
%             'NodeLabel', l, 'XData', scores(:, 1), 'YData', scores(:, 2), 'ZData', scores(:, 3));
    end
    title(sprintf('Cluster %d, Frequency = %.1f%%', icl, ...
        nMdlInCl(icl) / numel(allLmSt) * 100));
end
title(t, sprintf('%s model, graph of the limit sets', mdlName));

return

% Take a look at specific cluster

cloi = 3;  % Cluster of interest

fprintf('Unstable dimension for each fixed point in each model:\n\n')
us = {};
for i = 1:length(clIdx)
    if clIdx(i) == cloi
        idx1 = ceil(i / 2);
        idx2 = 1 + mod(i - 1, 2);
        curr = allLmSt{idx1, idx2};
        tmpIdx = find(vecnorm(curr.V, 2, 2) < 1e-4, 1);
        if tmpIdx
            origDim = curr.VUSDim(tmpIdx);
            tmp = ['Fixed points found: ' num2str(size(curr.V, 1)) ...
                '. Unstable dimension of the origin: ' num2str(origDim) ...
                '. All: ' num2str(sort(curr.VUSDim)') '. Model (' num2str([idx1, idx2]) ')'];
        else
            tmp = ['Fixed points found: ' num2str(size(curr.V, 1)) ...
                '. Origin: Not identified. All: ' num2str(sort(curr.VUSDim)') '. Model (' num2str([idx1, idx2]) ')'];
        end
        us = [us;{tmp}];
    end
end
us = sort(us);
disp(us)

% Plot limit set graphs for each model

nFigs = 0;
nx = 4; ny = 2;

for iFig = 1:nFigs
    f = figure();
    f.WindowState = 'maximized';
    t = tiledlayout(ny, nx);
    for y = 1:ny
        for x = 1:nx
            iSub = nx * (iFig - 1) + x;
            nexttile;
            curr = allLmSt{iSub, y};
            if iscell(curr)  % Limit Cycles
                nCy = numel(curr);
                if nCy == 1
                    viscircles([0 0], 1, 'Color', 'k');
                else
                    ang = linspace(0, 2 * pi, nCy + 1)';
                    ang = ang(1:end - 1);
                    viscircles(2 * [cos(ang) sin(ang)], 1, 'Color', 'k');
                end
                axis equal; xticks([]); yticks([]);
            else  % Fixed points
                adjMat = cellfun(@(x)~isempty(x), curr.E);
                c = distinguishable_colors(maxUSDim + 2);
                tmp = c(1, :); c(1, :) = c(2, :); c(2, :) = tmp;  % Use red for attractors
                cIdx = curr.VUSDim + 1;
                cIdx(cIdx > maxUSDim + 2) = maxUSDim + 2;
                c = c(cIdx, :);
                l = curr.VUSDim();
                
%                 % pick the nodes with connections
%                 idx = any(adjMat | adjMat');
%                 adjMat = adjMat(idx, idx);
%                 c = c(idx, :);
%                 l = l(idx);

                plot(digraph(adjMat), 'NodeColor', c, ...
                    'NodeLabel', l, 'Layout', 'force');
            end
            title(sprintf('Subject %02d, Model %d', iSub, y));
        end
    end
    title(t, sprintf('Model = %s, Graph of the limit sets', mdlName));
end

end