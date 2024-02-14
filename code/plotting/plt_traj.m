%% plt_traj - Plot the trajectories in PC space
%
% Ruiqi Chen, 2023/01/30
%
% Syntax:
%
% score = plt_traj(dat);
% score = plt_traj(dat, dim, dat_in_PC, nPC, coeff, mu);
% [score, coeff] = plt_traj(dat);
% [score, coeff, explained] = plt_traj(dat);  % Note |explained| is percentile
%
% Inputs:
%
% - |dat|:
%   - |dat_in_PC = false| (default): simulations, (nParcel, T, nSim) or (nParcel, T) matrix
%   - |dat_in_PC = true|: (T, nSim, K) matrix where |K| is number of PCs, scores in PC space
% - |dim|: 2 or 3, 2D or 3D visualization
% - |dat_in_PC|: whether |dat| is already in PC space, default is |false|
% - |nPC|: number of PCs to retain in the outputs, default is 10
% - |coeff|: if provided, |score| is computed by multiplying |dat| (after centering) with
%   |coeff| instead of calling |pca()|.
% - |mu|: if provided, |mu| is subtracted from |dat| instead of the mean of |dat|
function [score, varargout] = plt_traj(dat, dim, dat_in_PC, nPC, coeff, mu)

    if nargin < 6; mu = mean(dat, [2 3])'; end
    if nargin < 5; coeff = []; end
    if nargin < 4; nPC = 10; end
    if nargin < 3; dat_in_PC = false; end
    if nargin < 2; dim = 3; end
    assert(dim >= 2 && dim <= 3 && dim <= nPC)

    if ~dat_in_PC
        sims = permute(dat, [2 3 1]);
        [T, nSim, nParcel] = size(sims);
        sims = reshape(sims, [], nParcel);
        if ~isempty(coeff)
            coeff = coeff(:, 1:nPC);
            score = (sims - mu) * coeff;
            explained = [];
        else
            [coeff, score, ~, ~, explained, ~] = pca(sims, 'NumComponents', nPC);
        end
        score = reshape(score, T, nSim, []);
    else
        [T, nSim, K] = size(dat);
        coeff = eye(K, min([nPC K]));
        score = dat(:, :, 1:min([nPC K]));
        explained = [];
    end

    hold on;
    C = [0 0 1; 1 0 0];

    if dim == 2
        scatter(score(1, :, 1), score(1, :, 2), 24, ...
            C(1, :), 'filled');  % Starting points
        scatter(score(end, :, 1), score(end, :, 2), 80, ...
            C(2, :), 'filled');  % Ending points
        plot(score(:, :, 1), score(:, :, 2), 'k');
        xlabel("PC1"); ylabel("PC2");
    else
        scatter3(score(1, :, 1), score(1, :, 2), score(1, :, 3), 8, ...
            C(1, :), 'filled');  % Starting points
        scatter3(score(end, :, 1), score(end, :, 2), score(end, :, 3), 24, ...
            C(2, :), 'filled');  % Ending points
        plot3(score(:, :, 1), score(:, :, 2), score(:, :, 3), 'k');
        xlabel("PC1"); ylabel("PC2"); zlabel("PC3");
        view([-37.5, 30])
    end
    % legend(["Start", "End"])
    
    hold off;

    varargout = {};
    if nargout > 1
        varargout{1} = coeff;
    end
    if nargout > 2
        varargout{2} = explained;
    end

end