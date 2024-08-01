function [] = Censor_Models(mdlName, varargin)
%CENSOR_MODELS Plot the simulations of models that have been censored
%
%   Usage:
%     Censor_Models(mdlName, varargin)
%
%   Input:
%     mdlName: name of the model
%
%   Key-value inputs:
%     mdlFile: path to the model file
%     mtfFile: path to the motifs file. A boolean variable 'censored' will
%       be saved to this file.
%     censorFile: path to the censored models file. Should be a CSV with
%       columns 'Subject' (six-digit subjectID) and 'Session' (1 or 2).
%     figdir: path to the directory where the figures will be saved
%     vis: whether to save figures of simulations of censored models


%% Parameters

if nargin == 0 || isempty(mdlName)
    mdlName = 'HCP_Rest_FIX_Simple_Mdl200_sess';
end

p = inputParser;
addParameter(p, 'mdlFile', fullfile('data', [mdlName '.mat']), @ischar);
addParameter(p, 'mtfFile', fullfile('data', ['motifs_' mdlName '.mat']), @ischar);
addParameter(p, 'censorFile', fullfile('data', 'censored_models.csv'), @ischar);
addParameter(p, 'figdir', fullfile('figures', mdlName), @ischar);
addParameter(p, 'vis', true, @islogical);
parse(p, varargin{:});

mdlFile = p.Results.mdlFile;
mtfFile = p.Results.mtfFile;
censorFile = p.Results.censorFile;
figdir = p.Results.figdir;
vis = p.Results.vis;

if vis
    vfdir = fullfile(figdir, 'censored_models');
    if ~exist(vfdir, 'dir')
        mkdir(vfdir);
    end
end

% Visualization
nSimsShow = 120;  % Number of simulations to show
nSkipTR = 15;  % Number of initial TRs to skip
nRow = 2;  % Tile layout rows
nCol = 4;
nPlotsPerFig = nRow * nCol;
basefont = 14;  % Base font size


%% Generate censored data indicator

tbl = readtable(censorFile);
tbl.Subject = string(tbl.Subject);

load(mtfFile, 'sublist', 'allEq');
censored = false(size(allEq));
for i = 1:length(tbl.Subject)
    censored(sublist == tbl.Subject(i), tbl.Session(i)) = true;
end
save(mtfFile, 'censored', '-append');


%% Regenerate histogram for number & type of attractors

load(mtfFile, 'allLCMotif');
nEq = cellfun(@(x) size(x, 2), allEq);
nLC = cellfun(@(x) size(x, 2), allLCMotif);
nEqLC = string(nEq(:)) + "FP, " + string(nLC(:)) + "LC";

% Censoring
nEqLC_censored = nEqLC(censored(:));
nEqLC = nEqLC(~censored(:));

% Sort the categories based on popularity
nEqLC = categorical(nEqLC, "Ordinal", true);
[N, C] = histcounts(nEqLC);
[N, idx] = sort(N, 'descend');
C = C(idx);
P = N / numel(nEqLC);
nEqLC = reordercats(nEqLC, C);

fig = figure('Units', 'inches', 'Position', [0 0.5 10 5], 'Visible', 'off');
histogram(nEqLC, 'Normalization', 'probability');
hold on;
text(1:length(N), P, compose("%.3f", P), 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom');
xlabel('Number of stable fixed points (FP) and stable limit cycles (LC)');
ylabel('Proportion');
% title(sprintf('Distribution of attractors in %s (n = %d)', mdlName, numel(nEqLC)), ...
%     'Interpreter', 'none');
PrintAsSeen(fullfile(figdir, sprintf('nFPLC_%s', mdlName)), '-dpng', '-r300');
close(fig);


%% Visualize simulations of censored models

if ~vis
    return;
end

load(mdlFile, 'allMdl');
allMdl = allMdl(censored(:));
allEq = allEq(censored(:));
allLCMotif = allLCMotif(censored(:));

% Note: here models are ordered column-based but in the original CSV
% they were ordered row-based.
for i = 1:numel(allMdl)
    % Simulation
    disp(['Simulating model ' num2str(i) ' ...'])
    mdl = allMdl{i};
    inits = randn(size(mdl.Param{5}, 1), nSimsShow);
    [~, ~, sims] = GetEqLC(mdl, inits);

    % Create new figure
    if mod(i, nPlotsPerFig) == 1
        fig = figure('Units', 'inches', 'Position', [0 0.5 20 10.0208], ...
            'Visible', 'off');  % 27 inch screen maximized window size
        tlo = tiledlayout(fig, nRow, nCol, 'TileSpacing', 'compact', 'Padding', 'compact');
    end

    % Plot vector field
    nexttile();
    lgd = PlotAttractorLandscape(sims, allEq{i}, allLCMotif{i}, 3, nSkipTR);
    fontsize(gca, basefont, "points");
    title(sprintf('Identified: %s', nEqLC_censored(i)), ...
        'FontSize', basefont + 2);

    % Close old figure
    if i == numel(allMdl) || mod(i, nPlotsPerFig) == 0
        lgd.Layout.Tile = 'east';
        lgd.Visible = 'on';
        title(tlo, sprintf('Vector field & motifs of %s', mdlName), ...
            'Interpreter', 'none', 'FontSize', basefont + 4, 'FontWeight', 'bold');
        disp('Saving figure ...');
        PrintAsSeen(fullfile(vfdir, ...
            sprintf('Traj_censored_%s_%03d', mdlName, ceil(i / nPlotsPerFig))), ...
            '-dpng', '-r200');
        close(fig);
    end
end