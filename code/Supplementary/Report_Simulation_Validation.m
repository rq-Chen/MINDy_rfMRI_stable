function [mdlstr, nEqLC, Rsquared, Wcorr, Acorr, Dcorr, lmEig1] = Report_Simulation_Validation(mdlName, varargin)
%REPORT_SIMULATION_VALIDATION Reports about ground-truth simulation results
%
%   Usage:
%       [nEqLC, Rsquared, Wcorr, Acorr, Dcorr, lmEig1] = Report_Simulation_Validation(mdlName, varargin)
%
%   Positional inputs:
%     mdlName: name of the model (default: 'HCP_Rest_FIX_Simple_Mdl200_sess')
%
%   Key-value input:
%     mdlFile: name of the file containing the model (default: 'data/<mdlName>.mat')
%     datFile: name of the file containing the data (default: replace(mdlFile, '_Mdl', '_Deconv'))
%     derivType: type of derivatives (default: 'two-point')
%       - 'one-point': x'[t] = x[t + 1] - x[t];  (Default for MINDy_RC)
%       - 'two-point': x'[t] = (x[t + 2] - x[t]) / 2; (Default for MINDy_Simple)
%       - 'center': x'[t] = (x[t + 1] - x[t - 1]) / 2;
%     figdir: directory to save figures (default: 'figures/<mdlName>')
%     vis: whether to visualize the results for each model (default: true)
%     nSubs: number of subjects to process (default: 50). Note that we randomly pick one session.
%
%   Outputs:
%     mdlstr: names of the models, (1, 4) string
%     nEqLC: number of stable equilibria and limit cycles, (2, nSubs, 4) array
%     Rsquared: R squared of the fitting, (nSubs, 4) array
%     Wcorr: correlation between the W of all models and the first model, (nSubs, 4) array
%     Acorr: similar thing for the curvature A, (nSubs, 4) array
%     Dcorr: similar thing for the decay D, (nSubs, 4) array
%     lmEig1: largest magnitude of the eigenvalues of the linear model, (nSubs, 1) array
%
%   Ruiqi Chen, 09/27/2023
%
%   This script will compare these following models: MINDy fitted on data, "linear" MINDy fitted
%   on data, MINDy fitted on the noisy simulations of the obtained MINDy model, and MINDy fitted
%   on the noisy simulations of the "linear" MINDy model. The results will be ploted in a figure
%   and saved in the |figdir| directory. If |vis| is true, extra figures will be saved in the 
%   |validation| subdirectory with one figure for each subject.

init_prj('MINDy_Base_v1.0');


%% Parameters

if nargin == 0 || isempty(mdlName)
    mdlName = 'HCP_Rest_FIX_Simple_Mdl200_sess';
end

p = inputParser;
addParameter(p, 'mdlFile', fullfile('data', [mdlName '.mat']), @ischar);
addParameter(p, 'datFile', fullfile('data', [replace(mdlName, '_Mdl', '_Deconv') '.mat']), @ischar);
addParameter(p, 'mtfFile', fullfile('data', ['motifs_' mdlName '.mat']), @ischar);
addParameter(p, 'derivType', 'two-point', ...
    @(x) ismember(x, {'one-point', 'two-point', 'center'}));
addParameter(p, 'figdir', fullfile('figures', mdlName), @ischar);
addParameter(p, 'vis', true, @islogical);
addParameter(p, 'nSubs', 50, @isnumeric);
parse(p, varargin{:});

mdlFile = p.Results.mdlFile;
datFile = p.Results.datFile;
mtfFile = p.Results.mtfFile;
derivType = p.Results.derivType;
figdir = p.Results.figdir;
vis = p.Results.vis;
nSubs = p.Results.nSubs;

if ~exist(figdir, 'dir')
    mkdir(figdir);
end

if vis
    visdir = fullfile(figdir, 'validation');
    if ~exist(visdir, 'dir')
        mkdir(visdir);
    end
end

% Training parameters
nBatch = 5000;  % Default number of training batches
mindySimLen = 20;  % length of each segment of MINDy noisy simulation (make it short for stronger activation)
lmSimLen = 20;  % Similar thing for linear model


%% Select data

dF = matfile(datFile);
mF = matfile(mdlFile);

% Randomly select subjects
nTotalSubs = size(dF, 'allDat', 1);
nSubs = min(nSubs, nTotalSubs);
subIdx = randperm(nTotalSubs, nSubs);

load(mdlFile, 'runIdx');
if size(runIdx, 1) == 1
    runIdx = repmat(runIdx, nTotalSubs, 1);
end
nSess = size(runIdx, 2);

% Randomly select one session for each subject
sessIdx = randi(nSess, nSubs, 1);


%% Fit on surrogate data and linear MINDy

nEq = nan(nSubs, 3);  % number of stable equilibria (fitted on data/surrogate/linear)
nLC = nan(nSubs, 3);  % number of limit cycles (fitted on data/surrogate/linear)
R2 = nan(nSubs, 3);   % R squared of MINDy on surrogate, linear on data, MINDy on linear
load(mtfFile, 'allEq', 'allLCMotif');

for i = 1:nSubs
    disp("Processing subject " + i + "of" + nSubs + "...");
    nEq(i, 1) = size(allEq{subIdx(i), sessIdx(i)}, 2);
    nLC(i, 1) = size(allLCMotif{subIdx(i), sessIdx(i)}, 2);
    
    dat = dF.allDat(subIdx(i), runIdx{subIdx(i), sessIdx(i)});
    [~, ~, ~, matchedSpCov] = GetSurrogateData(dat);
    [~, eq, lc, R2(i, 1)] = GetStats(matchedSpCov, derivType, nBatch);
    nEq(i, 2) = size(eq, 2);
    nLC(i, 2) = size(lc, 2);
    
    % ---------Linear model-------- %

    % Get linear model
    [lmMdl, tmpX, tmpY] = MINDy_from_Deconv(dat, 'Linear', derivType, nBatch, false);
    tmpX = [tmpX{:}];
    tmpY = [tmpY{:}];

    % Get standard form x(t + 1) = Ax(t) and check for stability
    if ~isempty(lmMdl.Param{5})
        W=lmMdl.Param{5};
    else
        W=lmMdl.Param{1};
    end
    slope=lmMdl.Param{2};
    D=lmMdl.Param{6};
    A = W .* slope + diag(1 - D);  % Note: this is not the A for curvature
    nParcels = size(A, 1);
    lmEig1 = max(abs(eig(A)));

    % Get R squared and error magnitude for lmMdl
    lmErr = std(A * tmpX - tmpX - tmpY, 0, 2);
    R2(i, 2) = 1 - mean((lmErr .^ 2) ./ var(tmpY, 0, 2));

    % ---------MINDy on linear-------- %

    % Skip this if the linear model is unstable
    if lmEig1 > 1
        continue;
    end

    % Generate data
    nTR = size(tmpX, 2);
    nSegs = ceil(nTR / lmSimLen);
    tmpDat = cell(1, nSegs);
    for k = 1:nSegs
        tmp = nan(nParcels, lmSimLen);
        tmp(:, 1) = tmpX(:, randperm(nTR, 1));
        for l = 2:nTR
            tmp(:, l) = A * tmp(:, l - 1) + randn(nParcels, 1) .* lmErr;
        end
        tmpDat{k} = tmp;
    end
   
    [~, eq, lc, R2(i, 3)] = GetStats(tmpDat, derivType, nBatch);
    nEq(i, 3) = size(eq, 2);
    nLC(i, 3) = size(lc, 2);

end


%% Visualization

f = figure('Units', 'inches', 'Position', [1 2 10 6]);
t = tiledlayout(1, 2);

nexttile();
boxchart(R2);
xticklabels({'MINDy on surrogate data', 'linear MINDy on data', 'MINDy on linear MINDy'});
ylabel('R squared');
title('R squared of each model''s fit');

nexttile();
isMonoStable = (nEq == 1 & nLC == 0);
bar(sum(isMonoStable) ./ sum(~isnan(nEq)));
xticklabels({'MINDy on data', 'MINDy on surrogate data', 'MINDy on simulations of linear model'});
ylabel('Proportion');
title('Proportion of models with a single stable equilibrium');

fontsize(f, 14, 'points');
title(t, sprintf('Surrogate simulations (n = %d)', nSubs), 'FontSize', 16, 'FontWeight', 'bold');
PrintAsSeen(f, fullfile(figdir, ['surrogate_sims_' mdlName]), '-dpng', '-r300');


%% Simulation and fitting

% Models:
%  - MINDy fitted on data (mdl)
%  - MINDy fitted on the noisy simulations of the first model  (MonM)
%  - "linear" MINDy fitted on data (lmMdl)
%  - MINDy fitted on the noisy simulations of the second model (MonLM)
mdlstr = ["MINDy", "MINDy-on-MINDy", "Linear MINDy", "MINDy-on-linear-MINDy"];
nEqLC = nan(2, nSubs, 4);  % number of stable equilibria and limit cycles (page 3 not used)
Rsquared = nan(nSubs, 4);  % R squared of the fitting (towards either rfMRI or noisy simulations' derivative)
Wcorr = nan(nSubs, 4);  % correlation between the weight matrices of all models and the first model (page 1 not used)
Acorr = nan(nSubs, 4);  % similar thing for the curvature A (page 3 not used)
Dcorr = nan(nSubs, 4);  % similar thing for the decay D
lmEig1 = nan(nSubs, 1);  % largest magnitude of the eigenvalues of the linear model

for i = 1:nSubs
    
    disp("Processing subject " + i + "...");

    % Get data and model
    dat = dF.allDat(subIdx(i), runIdx{subIdx(i), sessIdx(i)});
    nParcels = size(dat{1}, 1);
    mdl = mF.allMdl(subIdx(i), sessIdx(i));
    mdl = mdl{1};

    % Get number of stable equilibria and limit cycles
    [mdlEq, mdlLC] = GetEqLC(mdl);
    nEqLC(1, i, 1) = size(mdlEq, 2);
    nEqLC(2, i, 1) = size(mdlLC, 2);

    % Get R squared and error magnitude for mdl
    switch derivType
        case 'two-point'
            tmpX = cellfun(@(x) x(:, 1:end - 2), dat, 'UniformOutput', false);
            tmpY = cellfun(@(x) (x(:, 3:end) - x(:, 1:end - 2)) / 2, dat, 'UniformOutput', false);
        case 'one-point'
            tmpX = cellfun(@(x) x(:, 1:end - 1), dat, 'UniformOutput', false);
            tmpY = cellfun(@(x) x(:, 2:end), dat, 'UniformOutput', false);
        case 'center'
            tmpX = cellfun(@(x) x(:, 2:end - 1), dat, 'UniformOutput', false);
            tmpY = cellfun(@(x) (x(:, 3:end) - x(:, 1:end - 2)) / 2, dat, 'UniformOutput', false);
    end
    tmpX = [tmpX{:}];
    tmpY = [tmpY{:}];
    pred = MINDyInt_00(mdl, tmpX, 1, 1, 0, 1);
    pred = squeeze(pred(:, 2, :) - pred(:, 1, :));
    mdlErr = std(tmpY - pred, 0, 2);
    Rsquared(i, 1) = 1 - mean((mdlErr .^ 2) ./ var(tmpY, 0, 2));

    % ---------MINDy on MINDy-------- %

    % Generate data
    nTR = size(tmpX, 2);
    nSegs = ceil(nTR / mindySimLen);
    tmpDat = cell(1, nSegs);
    for k = 1:nSegs
        tmpDat{k} = MINDyInt_00(mdl, tmpX(:, randperm(nTR, 1)), 1, 1, mdlErr, mindySimLen);
    end

    [MonM, MonMEq, MonMLC, Rsquared(i, 2), Wcorr(i, 2), Acorr(i, 2), Dcorr(i, 2)] = ...
        GetStats(tmpDat, derivType, nBatch, mdl);
    nEqLC(1, i, 2) = size(MonMEq, 2);
    nEqLC(2, i, 2) = size(MonMLC, 2);

    % ---------Linear model-------- %

    % Get linear model
    [lmMdl, tmpX, tmpY] = MINDy_from_Deconv(dat, 'Linear', derivType, nBatch, false);
    tmpX = [tmpX{:}];
    tmpY = [tmpY{:}];

    % Get standard form x(t + 1) = Ax(t) and check for stability
    if ~isempty(lmMdl.Param{5})
        W=lmMdl.Param{5};
    else
        W=lmMdl.Param{1};
    end
    slope=lmMdl.Param{2};
    D=lmMdl.Param{6};
    A = W .* slope + diag(1 - D);  % Note: this is not the A for curvature
    lmEig1(i) = max(abs(eig(A)));

    % Get R squared and error magnitude for lmMdl
    lmErr = std(A * tmpX - tmpX - tmpY, 0, 2);
    Rsquared(i, 3) = 1 - mean((lmErr .^ 2) ./ var(tmpY, 0, 2));

    % Correlation of parameters
    Wcorr(i, 3) = corr(W(:), mdl.Param{5}(:));
    Dcorr(i, 3) = corr(D(:), mdl.Param{6}(:));

    % ---------MINDy on linear-------- %

    % Skip this if the linear model is unstable
    if lmEig1(i) > 1
        continue;
    end

    % Generate data
    nTR = size(tmpX, 2);
    nSegs = ceil(nTR / lmSimLen);
    tmpDat = cell(1, nSegs);
    for k = 1:nSegs
        tmp = nan(nParcels, lmSimLen);
        tmp(:, 1) = tmpX(:, randperm(nTR, 1));
        for l = 2:nTR
            tmp(:, l) = A * tmp(:, l - 1) + randn(nParcels, 1) .* lmErr;
        end
        tmpDat{k} = tmp;
    end
   
    [MonLM, MonLMEq, MonLMLC, Rsquared(i, 4), Wcorr(i, 4), Acorr(i, 4), Dcorr(i, 4)] = ...
        GetStats(tmpDat, derivType, nBatch, mdl);
    nEqLC(1, i, 4) = size(MonLMEq, 2);
    nEqLC(2, i, 4) = size(MonLMLC, 2);
    
    % ---------Visualization-------- %

    if ~vis
        continue;
    end

    f = figure('Units', 'inches', 'Position', [0 0.5 20 10.0208], 'Visible', 'off');
    t = tiledlayout(2, 4);
    % MINDy W, MINDy motif1, MINDy motif2, lmMdl W;
    % MonM W, MonM motif1, MonM motif2, MonLM W

    % Compute most similar motifs between MINDy and MonM
    mdlLCMotif = LC2Motif(mdlLC);
    mdlMotifs = [mdlEq, mdlLCMotif{:}];
    MonMLCMotif = LC2Motif(MonMLC);
    MonMMotifs = [MonMEq, MonMLCMotif{:}];
    distMat = pdist2(mdlMotifs', MonMMotifs', "correlation");
    [mindist, idx] = min(distMat, [], 2);  % Plot motif 1, 2 for mdl and idx(1), idx(2) for MonM
    maxsim = 1 - mindist;

    nexttile;
    MyNetMatYeo(mdl.Param{5});
    title('MINDy connectivity W');
    subtitle(sprintf('R^2 = %.2f, %d FP, %d LC', Rsquared(i, 1), ...
        size(mdlEq, 2), size(mdlLC, 2)));
    
    for k = 1:2
        nexttile;
        if k > size(mdlMotifs, 2)  % Stable origin
            break
        end
        PlotYeoSurface(mdlMotifs(:, k));
        colorbar;
        title("MINDy motif " + k);
    end

    nexttile;
    MyNetMatYeo(W);
    title(sprintf('Linear MINDy W, offdiagcorr(W, MINDy W) = %.2f', Wcorr(i, 3)));
    subtitle(sprintf('R^2 = %.2f, max |\\lambda| = %.2f', Rsquared(i, 3), lmEig1(i)));

    nexttile;
    MyNetMatYeo(MonM.Param{5});
    title(sprintf('MINDy-on-MINDy W, offdiagcorr(W, MINDy W) = %.2f', Wcorr(i, 2)));
    subtitle(sprintf('R^2 = %.2f, %d FP, %d LC', Rsquared(i, 2), ...
        size(MonMEq, 2), size(MonMLC, 2)));
    
    for k = 1:2
        nexttile;
        if k > size(mdlMotifs, 2)  % Stable origin
            break
        end
        PlotYeoSurface(MonMMotifs(:, idx(k)));
        colorbar;
        title("MINDy-on-MINDy motif " + idx(k));
        subtitle(sprintf('Correlation with MINDy motif %d = %.2f', k, maxsim(k)));
    end

    nexttile;
    MyNetMatYeo(MonLM.Param{5});
    title(sprintf('MINDy-on-linear-MINDy W, offdiagcorr(W, MINDy W) = %.2f', Wcorr(i, 4)));
    subtitle(sprintf('R^2 = %.2f, %d FP, %d LC', Rsquared(i, 4), ...
        size(MonLMEq, 2), size(MonLMLC, 2)));
    
    title(t, sprintf("%s model ground truth simulations, subject %d, session %d", ...
        mdlName, subIdx(i), sessIdx(i)), 'Interpreter', 'none');
    subtitle(t, sprintf("MINDy simulation time = %d TRs, linear MINDy simulation time = %d TRs", ...
        mindySimLen, lmSimLen));
    PrintAsSeen(f, fullfile(visdir, sprintf('validation_%s_sub%d_sess%d', ...
        mdlName, subIdx(i), sessIdx(i))), '-dpng', '-r300');
    close(f);

end


%% Summary and visualization

stableOrig = squeeze(nEqLC(1, :, :) == 1 & nEqLC(2, :, :) == 0);
matchedDy = squeeze(all(nEqLC == nEqLC(:, :, 1)));  % (nSubs, 4)
matchType = categorical(stableOrig(:, [2 4]), [true false], {'Stable', 'Unstable'}) .* ...
    categorical(matchedDy(:, [2 4]), [true false], {'Matched', 'Mismatched'});
C = categorical(categories(matchType));

f = figure('Units', 'inches', 'Position', [0 0.5 20 10.0208]);
t = tiledlayout(1, 3);

nexttile;
N1 = histcounts(matchType(:, 1), C, 'Normalization', 'probability');
N2 = histcounts(matchType(:, 2), C, 'Normalization', 'probability');
bar(C, [N1' N2']);
ylabel('Proportion');
legend({'MINDy-on-MINDy', 'MINDy-on-linear-MINDy'});
title({'Stability of refitted models', 'and whether they match rfMRI MINDy dynamics'});

nexttile;
boxchart(Rsquared);
xticklabels(mdlstr);
ylabel('R squared');
title("R squared of each model's fit");

ax = nexttile;
x = {Wcorr(:, 2:end), Dcorr(:, 2:end), Acorr(:, 2:end)};
boxplotGroup(ax, x, 'primaryLabels', {'W', 'D', 'A'}, 'secondaryLabels', mdlstr(2:end), ...
    'Colors', lines(3), 'GroupType', 'betweenGroups', 'groupLabelType', 'Vertical');
% x = reshape([Wcorr(:, 2:end) Dcorr(:, 2:end) Acorr(:, 2:end)], [], 1);
% xgroups = repelem(["W" "D" "A"], nSubs, 3);
% xgroups = categorical(xgroups(:));
% cgroups = repmat(mdlstr(2:end), nSubs, 3);
% cgroups = categorical(cgroups(:));
% boxchart(xgroups, x, 'GroupByColor', cgroups, 'Notch', 'on');
% lgd = legend;
% lgd.Layout.Tile = 'east';
ylabel('Correlation');
title('Parameter correlation with rfMRI MINDy');

title(t, sprintf("%s model ground truth simulations (n = %d)", mdlName, nSubs), ...
    'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 16);
subtitle(t, sprintf("MINDy simulation time = %d TRs, linear MINDy simulation time = %d TRs", ...
    mindySimLen, lmSimLen));
PrintAsSeen(f, fullfile(figdir, sprintf('validation_%s', mdlName)), '-dpng', '-r300');


%% Clean up

% close all;
% rmpath('MINDy_Base_v1.0');

end


%% Helper functions

function [mdl, eq, lc, R2, Wcorr, Acorr, Dcorr] = GetStats(dat, derivType, nBatch, MINDyMdl)
%GETSTATS Get statistics of MINDy model fitted on |dat|

[mdl, tmpX, tmpY] = MINDy_from_Deconv(dat, 'Base', derivType, nBatch, false);
tmpX = [tmpX{:}];
tmpY = [tmpY{:}];

% Get number of stable equilibria and limit cycles
[eq, lc] = GetEqLC(mdl);

% Get R squared and error magnitude
pred = MINDyInt_00(mdl, tmpX, 1, 1, 0, 1);
pred = squeeze(pred(:, 2, :) - pred(:, 1, :));
err = std(tmpY - pred, 0, 2);
R2 = 1 - mean((err .^ 2) ./ var(tmpY, 0, 2));

% Correlation of parameters
if nargin > 3 && nargout > 4
    Wcorr = offdiagcorr(mdl.Param{5}, MINDyMdl.Param{5});
    Acorr = corr(mdl.Param{2}(:), MINDyMdl.Param{2}(:));
    Dcorr = corr(mdl.Param{6}(:), MINDyMdl.Param{6}(:));
end

end