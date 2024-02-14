function [eq, eqIdx, lc, lcIdx] = GetDataLimitSet(sims, EQ_TOL, RECUR_TOL, CL_TOL, LC_TOL, doPlot)
%GETDATALIMITSET Calculate the limit set according to simulated trajectories
%   Identify equilibria by scanning from front to back, looking for small
%   enough deviations. Identify limit cycles by computing the distance
%   between the states at each timepoint and the final state. Use the min
%   distance point to compute a putative period |T|, then check the
%   distance between different cycles.
%
%   Note that the algorithm requires at least one full cycle to identify
%   recurrence.
%
% Inputs:
%   - |sims|: (nParcel, T, nSims) matrix, output of |MINDyInt_00|
%   - |EQ_TOL|: tolerance for identifying equilibria (L-inf norm)
%   - |RECUR_TOL|: tolerance for identifying recurrence (L2 norm)
%   - |CL_TOL|: tolerance for grouping equilibria (L2 norm)
%   - |LC_TOL|: tolerance for grouping limit cycles (Mean L2 norm along the cycle)
%   - |doPlot|: Boolean, show plot for distance to the terminal point or not
%
% Outputs:
%   - |eq|: empty or (nParcels, ?) matrix, equilibria identified
%   - |eqIdx|: (nSims, 1) vector, the equilibria that each sim converges to
%       (0 if not converging to any eq)
%   - |lc|: empty or (1, ?) cell of (nParcels, ?) matrix, limit cycles
%   - |lcIdx|: (nSims, 1) vector, similar to |eqIdx|

% Inputs
if nargin < 6
    doPlot = false;
end
if nargin < 5 || isempty(LC_TOL)
    LC_TOL = 0.5;
end
if nargin < 4 || isempty(CL_TOL)
    CL_TOL = 1e-1;
end
if nargin < 3 || isempty(RECUR_TOL)
    RECUR_TOL = 0.5;
end
if nargin < 2 || isempty(EQ_TOL)
    EQ_TOL = 1e-6;
end

% Constants
[~, T, nSims] = size(sims);
% Identify convergence when |all(abs(dX(n)) < EQ_TOL)| for |n| in a range as least this long
EQ_TOL_LEN = min([T - 1 10]);
% Requiring at least this number of full cycles to identify recurrence
MIN_REC = 1;
MAX_STAYTIME_VAR_RATIO = 10;  % See below

% Identify point attractors
dx = sims(:, end - EQ_TOL_LEN + 1:end, :) - sims(:, end - EQ_TOL_LEN:end - 1, :);
reachEq = squeeze(all(abs(dx) < EQ_TOL, [1 2]));
if any(reachEq)
    % Final states of trajectories that converges to equilibria
    dup_eqs = squeeze(sims(:, end, reachEq));  % (nParcel, ?)

    eqIdx = zeros(nSims, 1);
    Idx = zeros(size(dup_eqs, 2), 1);

    % Add first equilibrium
    eq = dup_eqs(:, 1);
    Idx(1) = 1;

    % Scan through |dup_eqs| and add new equilibria
    for i = 2:size(dup_eqs, 2)
        [argval, argmin] = min(vecnorm(eq - dup_eqs(:, i)));
        if argval < CL_TOL
            Idx(i) = argmin;
        else
            eq = [eq dup_eqs(:, i)];
            Idx(i) = size(eq, 2);
        end
    end
    eqIdx(reachEq) = Idx;
else
    eq = [];
    eqIdx = zeros(nSims, 1);
end

if nargout <= 2
    return
end

% Identify recurrence
%
% Note that since dist2end(end, :) == 0, the last nonzero entry of
% "leaving" will be roughly when the state begin to leave x(T) in the last cycle,
% and the last nonzero entry of "approaching" is when it arrives around x(T)
% in the last cycle. Therefore, the time the state visits x(T) last time is
% about the center of the SECOND LAST nonzero entry of "approaching" and
% the LAST nonzero entry of "leaving".
dist2end = squeeze(vecnorm(sims - sims(:, end, :)));  % (T, nSims)
recurring = dist2end < RECUR_TOL;
approaching = (recurring(2:end, :) - recurring(1:end - 1, :) == 1);
leaving = (recurring(2:end, :) - recurring(1:end - 1, :) == -1);
hasLC = (sum(approaching) >= MIN_REC + 1) & (sum(leaving) >= MIN_REC);

% Exclude converging trajtories (usually caused by extremely slow spirals)
hasLC = hasLC & (~reachEq');

% Plot recurrence
if doPlot
    plot(dist2end);
    xlabel('Time (TR)');
    ylabel('Distance to final state (L2 norm)');
    title('Recurrence of the states');
    legend('Trajectories');
end

lcIdx = zeros(nSims, 1);
lc = {};

if any(hasLC)
    for i = 1:nSims
        if ~hasLC(i); continue; end
        
        a = find(approaching(:, i), 2, 'last');  % (second last, last)
        b = find(leaving(:, i), 1, 'last');
        
        % Exclude cases like:  --___------________________
        if (T - a(2) + 1) * 2 > MAX_STAYTIME_VAR_RATIO * (b - a(1) + 1)
            hasLC(i) = false;
            continue
        end

        a = a(1);

        [~, midx] = min(dist2end(a:b, i));
        currLC = sims(:, a + midx - 1:end, i);
        if isempty(lc)
            lc = {currLC};
            lcIdx(i) = 1;
        else
            for j = 1:length(lc)
                distMat = pdist2(currLC', lc{j}');  % (nT_curr, nT_lcj)
                if mean(min(distMat)) < LC_TOL  % mean distance of each lc{j} point to currLC
                    lcIdx(i) = j;
                    break;
                end
%                 period = min(size(lc{j}, 2), size(currLC, 2));
%                 a = lc{j}(:, 1:period); b = currLC(:, 1:period);
%                 distMat = pdist2(a', b');
%                 permDistMat = nan(period, period);
%                 for k = 1:period
%                     for l = 1:period  % column l: shifting the second loop by l
%                         permDistMat(k, l) = distMat(k, mod(k + l - 1, period) + 1);
%                     end
%                 end
%                 if min(mean(permDistMat)) < LC_TOL
%                     lcIdx(i) = j;
%                     break
%                 end
            end
            if lcIdx(i) == 0
                lc = [lc {currLC}];
                lcIdx(i) = length(lc);
            end
        end
    end
end

end