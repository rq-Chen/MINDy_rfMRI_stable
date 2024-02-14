%% Report_Param_Reliability.m - Analysis of parameters reliability
function [] = Report_Param_Reliability(mdlName, varargin)
%REPORT_RELIABILITY Reliability analysis of parameters
%
%   Usage:
%     Report_Param_Reliability(mdlName, ...)
%
%   Positional inputs:
%     mdlName: name of the model (default: 'HCP_Rest_FIX_Simple_Mdl200_sess')
%
%   Key-value inputs:
%     mdlFile: name of the file containing the model (dafault: 'data/<mdlName>.mat')
%     figdir: directory to save figures (default: 'figures/<mdlName>')


%% Parameters

if nargin == 0 || isempty(mdlName)
    mdlName = 'HCP_Rest_FIX_Simple_Mdl200_sess';
end

p = inputParser;
addParameter(p, 'mdlFile', fullfile('data', [mdlName '.mat']), @ischar);
addParameter(p, 'figdir', fullfile('figures', mdlName), @ischar);
parse(p, varargin{:});

mdlFile = p.Results.mdlFile;
figdir = p.Results.figdir;

if ~exist(figdir, 'dir')
    mkdir(figdir);
end

load(mdlFile, 'allMdl');
nSubs = size(allMdl, 1);
nSess = size(allMdl, 2);
assert(nSess > 1, 'Need at least two sessions for reliability analysis!');
if nSess > 2
    warning('More than two sessions, only using first two!');
end
basefont = 14;


%% Test-retest correlation across parameters (NOT across subjects!)

[withinW, acrossW] = GetParamCorr(cellfun(@(x) x.Param{5}, allMdl, 'UniformOutput', false));
[withinA, acrossA] = GetParamCorr(cellfun(@(x) x.Param{2}, allMdl, 'UniformOutput', false));
[withinD, acrossD] = GetParamCorr(cellfun(@(x) x.Param{6}, allMdl, 'UniformOutput', false));


%% ICC of the parameters across participants

ICCW = GetICC(cellfun(@(x) x.Param{5}, allMdl, 'UniformOutput', false));
ICCA = GetICC(cellfun(@(x) x.Param{2}, allMdl, 'UniformOutput', false));
ICCD = GetICC(cellfun(@(x) x.Param{6}, allMdl, 'UniformOutput', false));


%% Visualize parameter reliability

f = figure('Units', 'inches', 'Position', [0 0.5 12 10]);
t = tiledlayout(2, 2, 'Padding', 'compact');

nexttile();
h = boxplotGroup({[withinW withinD withinA], [acrossW acrossD acrossA]}, ...
    'PrimaryLabels', ["Within subject", "Across subject"], ...
    'SecondaryLabels', ["Connectivity", "Decay", "Curvature"], ...
    'Notch', 'on', ...
    'groupLabelType', 'Vertical', ...
    'Colors', [0 0 1; 0 0 0], 'GroupType', 'betweenGroups');
h.axis.XTickLabelRotation = 45; 
ylabel('Pearson correlation');
title('Correlation of parameters across sessions', 'FontSize', basefont + 2)

nexttile();
MyNetMatYeo(ICCW);
title('Intraclass correlation of connectivity', 'FontSize', basefont + 2);

nexttile();
PlotYeoSurface(ICCD);
colorbar;
title('Intraclass correlation of decay', 'FontSize', basefont + 2);

nexttile();
PlotYeoSurface(ICCA);
colorbar;
title('Intraclass correlation of curvature', 'FontSize', basefont + 2);

title(t, sprintf('Parameter reliability (n = %d)', nSubs), 'FontSize', basefont + 4, 'FontWeight', 'bold');
PrintAsSeen(f, fullfile(figdir, ['Param_Reliability_' mdlName]), '-dpng', '-r300');
close(f);

end


%% Utilities
function m = offdiagnanmean(M, varargin)
%OFFDIAGNANMEAN Mean of off-diagonal, non-NaN elements (if any; otherwise nan)
for i = 1:min(size(M))
    M(i, i) = nan;
end
m = mean(M, varargin{:}, "omitnan");
end


function ICCW = GetICC(W)
%GETRHO Get ICC of each weight across subjects
%
%   Inputs:
%     W: (nSubs, nSess) cell array of weights.
%
%   Outputs:
%     ICCW: same size as W{1}. ICC for each parameter.

sizeW = size(W{1});
nW = numel(W{1});
nSubs = size(W, 1);
% nSess = size(W, 2);
ICCW = nan(nW, 1);

allW = nan(nSubs, nW, 2);
for i = 1:2
    allW(:, :, i) = cell2mat(cellfun(@(x) x(:)', W(:, i), 'UniformOutput', false));
end

for i = 1:nW
    ICCW(i) = corr(allW(:, i, 1), allW(:, i, 2));
end
ICCW = reshape(ICCW, sizeW);

end


function [within, across, corrMat] = GetParamCorr(W)
%GETPARAMCORR Get correlation between weights in two sessions
%
%   Inputs:
%     W: (nSubs, nSess) cell array of weights.
%
%   Outputs:
%     within: (nSubs, 1) vector of within-subject correlation.
%     across: (nSubs, 1) vector of across-subject correlation,
%       averaged across sessions and other subjects.

nW = numel(W{1});
nSubs = size(W, 1);
% nSess = size(W, 2);

allW = nan(nW, nSubs, 2);
for i = 1:2
    allW(:, :, i) = cell2mat(cellfun(@(x) x(:), W(:, i)', 'UniformOutput', false));
end

corrMat = (corr(allW(:, :, 1), allW(:, :, 2)) + ...
    corr(allW(:, :, 2), allW(:, :, 1))) / 2;
within = diag(corrMat);
across = offdiagnanmean(corrMat)';

end