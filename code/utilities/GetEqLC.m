function [eq, lc, sims, eqIdx, lcIdx] = GetEqLC(mdl, inits, T)
%GETEQLC Get the stable fixed points and limit cycles of a model without input
%
%   Inputs:
%     - mdl: MINDy structure or EqFinderVec object
%     - inits: (nParcels, nSims) initial conditions (default 120 random inits)
%     - T: minimum simulation steps (default 1600). The function will keep
%         simulating until at least one attractor is found, up to 10T steps.
%
%   Outputs:
%     - eq: (nParcels, nEq) stable fixed points
%     - lc: (1, nLC) cell of (nParcels, tPeriod) limit cycles
%     - sims: (nParcels, nT * T + 1, nSims) simulations (nT >= 1)
%     - eqIdx: (nSims, 1) IDs of the eq that each sim converged to (0 if not
%         converging to any eq)
%     - lcIdx: (nSims, 1) IDs of the lc that each sim converged to (0 if not
%         converging to any lc)
%
%   The function will make sure the number of stable fixed points are even
%   by computing a correlation matrix and adding the missing one in a pair.

    % Handle inputs
    if isstruct(mdl)
        f = @(m, x0, t1) MINDyInt_00(m, x0, 1, 1, 0, t1);
        nParcels = length(mdl.Param{5});
    elseif isa(mdl, 'EqFinderVec')
        f = @(m, x0, t1) m.run(x0, t1);
        nParcels = mdl.nParcels;
    else
        error('Input |mdl| must be a MINDy structure or EqFinderVec object')
    end

    if nargin < 3
        T = 1600;
    end
    if nargin < 2
        inits = randn(nParcels, 120);
    end
    cosThres = 0.95;  % Minimum cosine similarity to be considered as a pair
    TreatedAsOrigin = 1e-6;  % Treat a fixed point as the origin if mean square < this

    % Simulate |T| steps each time until one attractor was found
    for nT = 1:10
        currsims = f(mdl, inits, T);
        if nT == 1
            sims = currsims;  % (nParcels, T + 1, nSims)
        else
            sims = cat(2, sims, currsims(:, 2:end, :));
        end
        [eq, eqIdx, lc, lcIdx] = GetDataLimitSet(sims);
        nAtt = size(eq, 2) + numel(lc);
        if nAtt == 0  % Need more time
            inits = squeeze(sims(:, end, :));
        else
            break
        end
    end
    
    % Add missing fixed points (NOT suitable for asymmetric models)

    % Check the origin first
    idx = find(mean(eq .^ 2) < TreatedAsOrigin);
    if ~isempty(idx)
        origin_stable = true;
        eq(:, idx) = [];
        warning('Stable origin identified');
    else
        origin_stable = false;
    end
    % Add symmetric fixed points
    if mod(size(eq, 2), 2)
        distMat = squareform(pdist(eq', 'cosine'));
        distMat = abs(1 - distMat);  % 0 orthogonal, 1 the same or reflected
        distMat = distMat - diag(diag(distMat));  % Set diagonal to 0
        tmp = all(distMat < cosThres);
        if sum(tmp)
            eq = [eq -eq(:, tmp)];  % Adding from the back thus |eqIdx| is not influenced
            warning('Adding %d symmetric fixed points', sum(tmp))
        end
    end
    % Add the origin back
    if origin_stable
        eq = [zeros(size(eq, 1), 1) eq];
    end
end

