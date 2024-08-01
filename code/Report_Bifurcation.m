function [nEqLCw, nEqLCb, nEqLCOrig] = Report_Bifurcation(mdlName, varargin)
%REPORT_BIFURCATION Analysis of convex combination of two models
%
%   Usage:
%     [nEqLCw, nEqLCb, nEqLCOrig] = Report_Bifurcation(mdlName, varargin);
%
%   Positional input:
%     mdlName: name of the model to be analyzed
%
%   Key-value inputs:
%     mdlFile: name of the file containing the model (default: 'data/<mdlName>.mat')
%     mtfFile: name of the file containing the motifs (default: 'data/motifs_<mdlName>.mat')
%     figdir: name of the directory where to save the figures (default: 'figures/<mdlName>')
%     subID: ID of the subject to analyze in details (set as "" to skip)
%     nSamp: number of linearly combined models to use for each distribution (default: 500)
%
%   Outputs:
%     nEqLCw: number of fixed points and limit cycles for linearly combined models within subject
%     nEqLCb: number of fixed points and limit cycles for linearly combined models between subjects
%     nEqLCOrig: number of fixed points and limit cycles for the original models

if nargin == 0 || isempty(mdlName)
    mdlName = 'HCP_Rest_FIX_Simple_Mdl200_sess';
end

p = inputParser;
addParameter(p, 'mdlFile', fullfile('data', [mdlName '.mat']), @ischar);
addParameter(p, 'mtfFile', fullfile('data', ['motifs_' mdlName '.mat']), @ischar);
addParameter(p, 'figdir', fullfile('figures', mdlName), @ischar);
addParameter(p, 'subID', "210011", @(x) ischar(x) || isstring(x));
addParameter(p, 'nSamp', 500, @isnumeric)
parse(p, varargin{:});

mdlFile = p.Results.mdlFile;
mtfFile = p.Results.mtfFile;
figdir = p.Results.figdir;
subID = string(p.Results.subID);
nSamp = p.Results.nSamp;

if ~exist(figdir, 'dir')
    mkdir(figdir);
end

% Other parameters
nCmb = 6;  % Number of linear combinations for each pair of models
w = linspace(1, 0, nCmb);  % Weights (w * mdl1 + (1 - w) * mdl2)
nrows = 2;  % Number of rows for the figure
nAllcols = 4;  % The first column is used for the histogram
nfcols = nAllcols - 1;
assert(nrows * nfcols >= nCmb, 'Not enough panels for selected number of combinations');
basefont = 14;

mF = matfile(mdlFile);
load(mdlFile, 'sublist');
sublist = extractBetween(sublist, "sub", "Y");


%% Detailed analysis of a subject

if ~strcmp(subID, "")

    % Load data
    iSub = find(strcmp(sublist, subID));
    if isempty(iSub)
        error('Subject not found in the input file');
    end
    mdls = mF.allMdl(iSub, :);

    % Plot
    figure('Units', 'inches', 'Position', [0 0.5 20 10.0208]);
    t = tiledlayout(nrows, nAllcols, 'TileSpacing', 'compact', 'Padding', 'compact');
    for i = 1:nCmb
        nexttile(i + ceil(i / nfcols));  % Right nfcols columns
%         nexttile(tilenum(t, ceil(i / nfcols), 1 + mod(i - 1, nfcols)));  % Left nfcols columns
        mdl = EqFinderVec(mdls, [w(i) 1 - w(i)]);
        [eq, lc, sims] = GetEqLC(mdl);
        lcmtf = LC2Motif(lc);
        if i == 1
            [~, coeff, ~, mu] = myPCA(sims);
        end
        lgd = PlotAttractorLandscape(sims, eq, lcmtf, 3, [], coeff, mu);
        view(-64, 29);
        fontsize(gca, basefont, "points");
        title(sprintf('%.2f model 1 + %.2f model 2', w(i), 1 - w(i)), 'FontSize', basefont + 2);
    end
    lgd.Layout.Tile = 'east';
    lgd.Visible = 'on';
    lgd.AutoUpdate = 'off';
    title(t, ['Linear interpolations of ' mdlName ' models'], 'Interpreter', 'none', ...
        'FontSize', basefont + 4, 'FontWeight', 'bold');
    subtitle(t, "Example subject: " + subID, 'FontSize', basefont + 2, 'FontAngle', 'italic');

end


%% Taxonomy of attractors for linearly combined models

load(mdlFile, 'allMdl');

% Remove all subjects with at least one model censored
try
    load(mtfFile, 'censored');
    allMdl = allMdl(~any(censored, 2), :);
catch
    warning('Did not find censored data index. Assuming no censored data.')
end
nSubs = size(allMdl, 1);
nMdls = numel(allMdl);

% Linearly combined models
nEq = zeros(nSamp, 2);  % within, between
nLC = zeros(nSamp, 2);
for i = 1:nSamp
    disp("Generating sample " + i + " of " + nSamp + " ...")
    wi = rand();
    % Within subject
    mdl1 = EqFinderVec(allMdl(randi(nSubs), :), [wi 1 - wi]);
    % Across subjects
    mdl2 = EqFinderVec(allMdl(randperm(nMdls, 2)), [wi 1 - wi]);
    mdls = {mdl1, mdl2};
    for j = 1:2
        [eq, lc] = GetEqLC(mdls{j});
        nEq(i, j) = size(eq, 2);
        nLC(i, j) = size(lc, 2);
    end
end
nEqLC = string(nEq) + "FP, " + string(nLC) + "LC";
nEqLCw = categorical(nEqLC(:, 1));
nEqLCb = categorical(nEqLC(:, 2));

% Original models
load(mtfFile, 'allEq', 'allLCMotif');
if exist("censored", "var")
    allEq = allEq(~any(censored, 2), :);
    allLCMotif = allLCMotif(~any(censored, 2), :);
end

nEqOrig = cellfun(@(x) size(x, 2), allEq);
nLCOrig = cellfun(@(x) size(x, 2), allLCMotif);
nEqLCOrig = string(nEqOrig(:)) + "FP, " + string(nLCOrig(:)) + "LC";
nEqLCOrig = categorical(nEqLCOrig);


%% Visualize the taxonomy

% Rename categories with less than |pthres| frequency to "others"
pthres = 5e-3;
[N, C] = histcounts(nEqLCOrig, 'Normalization', 'probability');
[~, idx] = sort(N, 'descend');
N = N(idx); C = C(idx);
Ckeep = C(N >= pthres);
Call = unique([categories(nEqLCOrig); categories(nEqLCw); categories(nEqLCb)]);
Crmv = Call(~ismember(Call, Ckeep));
rmvName = sprintf("Others (p < %.0g)", pthres);
C = [Ckeep'; rmvName];
C = categorical(C, C, C);
C = reordercats(C, string(C));
nEqLCOrig = mergecats(nEqLCOrig, Crmv, rmvName);
nEqLCw = mergecats(nEqLCw, Crmv, rmvName);
nEqLCb = mergecats(nEqLCb, Crmv, rmvName);

% Count
N1 = histcounts(nEqLCOrig, C, 'Normalization', 'probability');
N2 = histcounts(nEqLCw, C, 'Normalization', 'probability');
N3 = histcounts(nEqLCb, C, 'Normalization', 'probability');
N = [N1; N2; N3]';

% Horizontal bar plot
if subID ~= ""
    nexttile(1, [nrows 1]);
    % nexttile(nAllcols, [nrows, 1]);  % On the right
else
    figure;
end
hold on;
b = barh(C, N);
for i = 1:numel(b)
    text(b(i).YEndPoints, b(i).XEndPoints, compose(" %.3g", b(i).YData), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end
xlim([8e-4 1]);
ylabel('Number of stable fixed points (FP) and stable limit cycles (LC)');
xlabel('Frequency');
legend(["Original" "Interpolated within participant" "Interpolated across participants"], ...
    "Location", "best");
fontsize(gca, basefont, "points");
title('Distribution of the types of attractors', 'Interpreter', 'none', ...
    'FontSize', basefont + 2);
subtitle(sprintf('observed (%d models), interpolated (%d models)', nMdls, nSamp), ...
    "FontSize", basefont, 'FontAngle', 'italic');

% Compatibility with the example plots
if subID ~= ""
    lgd = gcf().Children.Children(4);
    lgd.Layout.Tile = 'east';
    lgd.Visible = 'on';
end

% Save figure
PrintAsSeen(fullfile(figdir, ['bifurcations_' mdlName]), '-dpng', '-r300');
% PrintAsSeen(fullfile(figdir, ['bifurcations_' mdlName]), '-dpdf', '-vector');


end