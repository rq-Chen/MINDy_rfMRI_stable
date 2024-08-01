function [scaling, simDat, simFC, datFC] = GetOptimalNoise(allMdl, allDat, runIdx, nInits, mins, maxs, initx0, option)
%GETOPTIMALNOISE Search for an optimal scaling of the noise to match the FC
%
% [scaling, simDat, simFC, datFC] = GetOptimalNoise(allMdl, allDat, runIdx, ...
%     nInits, mins, maxs, initx0, option)
%
% If `mins` equals `maxs`, the function will simulate data at the given scaling,
% otherwise it will search for the optimal scaling in the range [`mins`, `maxs`]
% and return the simulated data at the optimal scaling.
%
% If `initx0` is `true`, the model simulation will be initialized with the first
% timepoint of the data in each run in the order as used in training, otherwise
% it uses a random timepoint.
%
% We compute a single scaling factor which maximally matched the empirical and
% simulated FC across all subjects and sessions.
%
% To check the optimization process, set `option = optimset('Display', 'iter', 
% 'PlotFcns', {@optimplotx, @optimplotfval}, ...)`.

if nargin < 4 || isempty(nInits)
    nInits = 2;  % Number of runs to simulate per model
end
if nargin < 5 || isempty(mins)
    mins = 0.1;
end
if nargin < 6 || isempty(maxs)
    maxs = 1;
end
if nargin < 7 || isempty(initx0)
    initx0 = false;  % Initialize with a random timepoint
end
if nargin < 8 || isempty(option)
    option = optimset('Display', 'notify', 'TolX', 1e-3, 'MaxFunEvals', 50);
end

[nSubs, nSess] = size(allMdl);

% Compute data FC
datFC = cell(nSubs, nSess);
for i = 1:nSubs
    for j = 1:nSess
        datFC{i, j} = corrcoef([allDat{i, runIdx{i, j}}]');
    end
end
mDatFC = mean(cat(3, datFC{:}), 3);

% Search for the optimal scaling
simDat = cell(nSubs, nSess, nInits);
simFC = cell(nSubs, nSess);
if mins == maxs  % Just simulate data at the given scaling
    scaling = mins;
    GetFCErr(scaling);
else
    scaling = fminbnd(@GetFCErr, mins, maxs, option);
end

% Compute the mismatch between the model FC and the data FC
function FCErr = GetFCErr(s)
    for i1 = 1:nSubs
        for j1 = 1:nSess
            nRuns = numel(runIdx{i1, j1});
            for k1 = 1:nInits
                if initx0
                    idx = runIdx{i1, j1}(mod(k1 - 1, nRuns) + 1);
                    nTp = size(allDat{i1, idx}, 2);
                    x0 = allDat{i1, idx}(:, 1);
                else
                    idx = runIdx{i1, j1}(randi(nRuns));
                    nTp = size(allDat{i1, idx}, 2);
                    x0 = allDat{i1, idx}(:, randi(nTp));
                end
                simDat{i1, j1, k1} = MINDyInt_00(allMdl{i1, j1}, x0, 1, 1, ...
                    allMdl{i1, j1}.RMSE * s, nTp - 1);
            end
            simFC{i1, j1} = corrcoef([simDat{i1, j1, :}]');
        end
    end
    mSimFC = mean(cat(3, simFC{:}), 3);
    FCErr = mean((trilv(mSimFC) - trilv(mDatFC)) .^ 2);
end

end


%% Utility function

function [A] = trilv(A)
    A = squareform(tril(A, -1))';
end