%% Report_Attractor_Reliability.m - Analysis of parameters & motif reliability
function [] = Report_Attractor_Reliability(mdlName, varargin)
%REPORT_ATTRACTOR_RELIABILITY Reliability analysis of attractors
%
%   Usage:
%     Report_Reliability(mdlName, ...)
%
%   Positional inputs:
%     mdlName: name of the model (default: 'HCP_Rest_FIX_Simple_Mdl200_sess')
%
%   Key-value inputs:
%     mdlFile: name of the file containing the model (dafault: 'data/<mdlName>.mat')
%     mtfFile: name of the file containing the motifs (default: 'data/motifs_<mdlName>.mat')
%     datFile: name of the file containing the data (default: replace(mdlFile, '_Mdl', '_Deconv'))
%     figdir: directory to save figures (default: 'figures/<mdlName>')
%
%   Currently |datFile| is not used and we only compute motif reliability.


%% Parameters

if nargin == 0 || isempty(mdlName)
    mdlName = 'HCP_Rest_FIX_Simple_Mdl200_sess';
end

p = inputParser;
addParameter(p, 'mdlFile', fullfile('data', [mdlName '.mat']), @ischar);
addParameter(p, 'mtfFile', fullfile('data', ['motifs_' mdlName '.mat']), @ischar);
addParameter(p, 'figdir', fullfile('figures', mdlName), @ischar);
addParameter(p, 'LCGhostOnly', true, @islogical);
addParameter(p, 'exclude_n_LC', 0, @isnumeric);  % Exclude models with more than this number of LCs
parse(p, varargin{:});

mdlFile = p.Results.mdlFile;
mtfFile = p.Results.mtfFile;
figdir = p.Results.figdir;
LCGhostOnly = p.Results.LCGhostOnly;
exclude_n_LC = p.Results.exclude_n_LC;

if ~exist(figdir, 'dir')
    mkdir(figdir);
end

basefont = 14;


%% Load data

load(mtfFile, 'allEq', 'allLCMotif');
try  % For compatibility with outdated data files
    size(allLCMotif);
catch
    load(mtfFile, 'allLC');
    allLCMotif = cellfun(@LC2Motif, allLC, 'UniformOutput', false);
end
% Load model
load(mdlFile, 'allMdl');

% Exclude censored data
try
    load(mtfFile, 'censored');
    allEq = allEq(~any(censored, 2), :);
    allLCMotif = allLCMotif(~any(censored, 2), :);
    allMdl = allMdl(~any(censored, 2), :);
catch
    warning('Did not find censored data index. Assuming no censored data.')
end

if LCGhostOnly  % Reduce the bias of LC having too many motifs
    for i = 1:numel(allLCMotif)
        if size(allLCMotif{i}, 2) == 1
            allLCMotif{i}{1} = allLCMotif{i}{1}(:, [1 3]);
        elseif size(allLCMotif{i}, 2) > 1
            allLCMotif{i} = cellfun(@(x) x(:, 1), allLCMotif{i}, 'UniformOutput', false);
        end
    end
end

% Remove models with too many LCs
nLC = cellfun(@(x) size(x, 2), allLCMotif);
if exclude_n_LC
    keepSub = all(nLC <= exclude_n_LC, 2);
else
    keepSub = true(size(nLC, 1), 1);
end
allEq = allEq(keepSub, :);
allLCMotif = allLCMotif(keepSub, :);
allMdl = allMdl(keepSub, :);
[nSubs, nSess] = size(allEq);


%% Attractor consistency

% Attractor types and consistency across sessions
nEq = cellfun(@(x) size(x, 2), allEq);
nLC = cellfun(@(x) size(x, 2), allLCMotif);
nEqLC = string(nEq) + "FP, " + string(nLC) + "LC";
nEqLC = categorical(nEqLC);

% Remove trivial attractors
allEq = cellfun(@(x) x(:, any(x)), allEq, 'UniformOutput', false);

% Combine LC cells
allLCMotif = cellfun(@(x) [x{:}], allLCMotif, 'UniformOutput', false);

% Combine all motifs
allV = cellfun(@(x, y) [x y], allEq, allLCMotif, "UniformOutput", false);

% Get maximum correlation between motifs across all subjects and sessions
MaxCorr = nan(nSubs, nSubs, nSess, nSess);
for i = 1:nSubs
    for j = 1:nSubs
        for k = 1:nSess
            for l = 1:nSess
                if k == l  % Only compare across sessions
                    continue
                end
                if isempty(allV{i, k}) || isempty(allV{j, l})  % Trivial attractor only
                    continue
                elseif i > j || (i == j && k > l)
                    continue
                    % MaxCorr(i, j, k, l) = MaxCorr(j, i, l, k);
                elseif i == j && k == l
                    continue
                    % MaxCorr(i, j, k, l) = 1;
                else
                    % Note: model is symmetric so we take the absolute value
                    MaxCorr(i, j, k, l) = ...
                        max(abs(corr(allV{i, k}, allV{j, l})), [], 'all');
                end
            end
        end
    end
end

% % Fisher's transformation
% MaxCorr = atanh(MaxCorr);

WithinOrAcross = false(nSubs, nSubs, nSess, nSess);
for i = 1:nSubs
    WithinOrAcross(i, i, :, :) = true;
end
WithinOrAcross = categorical(WithinOrAcross, [true, false], {'Same person', 'Different person'});
AttTypeConsistency = false(nSubs, nSubs, nSess, nSess);
for i = 1:nSubs
    for j = 1:nSubs
        for k = 1:nSess
            for l = 1:nSess
                AttTypeConsistency(i, j, k, l) = ...
                    nEqLC(i, k) == nEqLC(j, l);
            end
        end
    end
end
AttTypeConsistency = categorical(AttTypeConsistency, [true, false], {'Same dynamics', 'Different dynamics'});
MaxCorr = MaxCorr(:);
WithinOrAcross = WithinOrAcross(:);
AttTypeConsistency = AttTypeConsistency(:);


%% Visualize motif reliability

f = figure('Units', 'inches', 'Position', [0 0.5 14 6]);
t = tiledlayout(1, 2);

nexttile(2);
aov = anova(table(MaxCorr, WithinOrAcross, AttTypeConsistency), 'MaxCorr', ...
    'ModelSpecification', 'interactions');
tmp = stats(aov);
boxchart(aov, ["WithinOrAcross", "AttTypeConsistency"]);

ylabel('Correlation');
% ylabel("Fisher's Z_{r}")
lgd = legend;
lgd.AutoUpdate = 'off';

% if tmp{"WithinOrAcross", "pValue"} < 0.05
%     hold on;
%     plot([1 2], [0.1 0.1], '-', 'Color', [0.5 0.5 0.5]);
%     % plot([1 1], [0.1 0.11], '-k');
%     % plot([2 2], [0.1 0.11], '-k');
%     text(1.4, 0.11, '* p < 0.05', 'Color', [0.5 0.5 0.5], 'VerticalAlignment', 'baseline');
% end

fontsize(gca, basefont, "points");
title({'Dominant Attractor Similarity (DAS)', ''});
% title(sprintf('MaxCorr ~ 1 + WithinOrAcross * AttTypeConsistency (n = %d)', nSubs), ...
%     'FontSize', basefont + 2);
% subtitle(sprintf('P(WithinOrAcross) = %.2g, P(AttTypeConsistency) = %.2g, P(interaction) = %.2g', ...
%     tmp{"WithinOrAcross", "pValue"}, tmp{"AttTypeConsistency", "pValue"}, ...
%     tmp{"WithinOrAcross:AttTypeConsistency", "pValue"}), "FontSize", basefont, "FontAngle", "italic");

% title(t, 'Maximum inter-session correlation between motifs', 'FontSize', basefont + 4, 'FontWeight', 'bold');
% subtitle(t, mdlName, 'Interpreter', 'none', "FontSize", basefont + 2, "FontAngle", "italic");
% PrintAsSeen(f, fullfile(figdir, ['Motif_Reliability_' mdlName]), '-dpng', '-r300');
% close(f);


%% Limit cycle speed distribution

% Remove the added ghost for 1LC models
if LCGhostOnly
    for i = 1:numel(allLCMotif)
        if nLC(i) == 1
            allLCMotif{i} = allLCMotif{i}(:, 1);
        end
    end
end

spRange = cellfun(@GetSpeedRange, allMdl(nLC > 0), allLCMotif(nLC > 0), "UniformOutput", false);
spRange = [spRange{:}];

nexttile(1);
histogram(log10(spRange), 'Normalization', 'probability');
hold on;
yl = ylim();
tmp = median(spRange);
xline(log10(tmp), '--k');
text(log10(tmp), yl(1) * 0.2 + yl(2) * 0.8, sprintf('  median = %.1f', tmp), ...
    "FontSize", basefont);
tmp = [1 2 4 10 20 40 100 200];
xticks(log10(tmp));
xticklabels(tmp);
xlabel('Ratio between max and min speed');
ylabel('Frequency');
fontsize(gca(), basefont, 'points');
title({sprintf('Distribution of speed range on limit cycles (n = %d)', numel(spRange)), ''}, 'FontSize', basefont + 2, ...
    'FontWeight', 'bold');

lgd.Layout.Tile = 'east';
set(findall(gcf,'type','text'), 'FontName','helvetica');
PrintAsSeen(f, fullfile(figdir, ['Motif_Reliability_' mdlName]), '-dsvg', '-vector');
PrintAsSeen(f, fullfile(figdir, ['Motif_Reliability_' mdlName]), '-dpdf', '-vector');
save(fullfile('data', ['LCspRange_' mdlName '.mat']), 'spRange');

end


%% Utilities
function m = offdiagnanmean(M, varargin)
%OFFDIAGNANMEAN Mean of off-diagonal, non-NaN elements (if any; otherwise nan)
for i = 1:min(size(M))
    M(i, i) = nan;
end
m = mean(M, varargin{:}, "omitnan");
end

function r = GetSpeedRange(mdl, init)
    [~, lc] = GetEqLC(mdl, init);
    r = nan(size(lc));
    for i = 1:numel(lc)
        v = squeeze(vecnorm(lc{i}(:, 2:end) - lc{i}(:, 1:end-1)));
        r(i) = max(v) / min(v);
    end
end