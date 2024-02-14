function [sublist, allEq, allLCMotif] = Report_Attractor_Landscape(mdlName, varargin)
%REPORT_ATTRACTOR_LANDSCAPE Visualize vector field and attractor motifs of the model.
%
%   Usage:
%       [sublist, allEq, allLCMotif] = REPORT_ATTRACTOR_LANDSCAPE(mdlName, varargin)
%
%   Positional input:
%     mdlName: name of the model (default: 'HCP_Rest_FIX_Simple_Mdl200_sess')
%
%   Key-value input:
%     mdlFile: name of the file containing the model (default: 'data/<mdlName>.mat')
%     outFile: name of the file to save the extracted motifs (default: 'data/motifs_<mdlName>.mat')
%     figdir: directory to save the figures (default: 'figures/<mdlName>')
%     vis: whether to visualize the vector field and motifs (default: true)
%     nSubs: number of subjects to process (default: 0, all subjects)
%
% Ruiqi Chen, 09/02/2023
%
% This script visualizes the vector field and attractor motifs of the model.
% Figures are saved to |figdir| and the extracted motifs are saved to |outFile|.
%
% |outFile| (and the function output) contains the following variables:
%   - |sublist|: subject IDs. (nSubs, 1) string.
%   - |allEq|: stable fixed points. (nSubs, nSess) cell array of (nParcels, nEq)
%     matrices. All found stable fixed points are stored, including reflections.
%   - |allLCMotif|: attractor motifs. (nSubs, nSess) cell array of (1, nLC) cell.
%     Each cell contains a (nParcels, 4) matrix, where each column is:
%       1. "Positive semi-major axis extreme";
%       2. "Positive semi-minor axis extreme";
%       3. "Negative semi-major axis extreme";
%       4. "Negative semi-minor axis extreme".
%
% The limit cycle motifs are extracted using the default method of |LC2Motif()|.


%% Parameters

if nargin == 0 || isempty(mdlName)
    mdlName = 'HCP_Rest_FIX_Simple_Mdl200_sess';
end

p = inputParser;
addParameter(p, 'mdlFile', fullfile('data', [mdlName '.mat']), @ischar);
addParameter(p, 'outFile', fullfile('data', ['motifs_' mdlName '.mat']), @ischar);
addParameter(p, 'figdir', fullfile('figures', mdlName), @ischar);
addParameter(p, 'vis', true, @islogical);
addParameter(p, 'nSubs', 0, @isnumeric);
parse(p, varargin{:});

mdlFile = p.Results.mdlFile;
outFile = p.Results.outFile;
figdir = p.Results.figdir;
vis = p.Results.vis;
nSubs = p.Results.nSubs;

if ~exist(figdir, 'dir')
    mkdir(figdir);
end
vfdir = fullfile(figdir, 'vector_field');
if ~exist(vfdir, 'dir')
    mkdir(vfdir);
end

% Visualization
nSimsShow = 60;  % Number of simulations to show
nSkipTR = 15;  % Number of initial TRs to skip
nRow = 2;  % Tile layout rows
nCol = 4;
nPlotsPerFig = nRow * nCol;
basefont = 14;  % Base font size


%% Initialization

load(mdlFile, 'allMdl');
if nSubs == 0
    nSubs = size(allMdl, 1);
else
    allMdl = allMdl(1:nSubs, :);
end
try
    load(mdlFile, 'sublist');
    sublist = extractBetween(sublist, "sub", "Y");
    sublist = sublist(1:nSubs);
catch
    sublist = "#" + string((1:nSubs)');
end

nSess = size(allMdl, 2);
allEq = cell(nSubs, nSess);
allLCMotif = cell(nSubs, nSess);


%% Show vector fields and motifs

nPlots = 0;
for iSub = 1:nSubs
    for iSess = 1:nSess

        fprintf('Processing subject %s, model %d ...\n', sublist(iSub), iSess);
        mdl = allMdl{iSub, iSess};

        % Get trajectories and attractors
        [eq, lc, sims] = GetEqLC(mdl);

        % Get motifs
        lcMotif = LC2Motif(lc);

        % Save motifs
        allEq{iSub, iSess} = eq;
        allLCMotif{iSub, iSess} = lcMotif;

        if ~vis
            continue;
        end

        % Create new figure
        if mod(nPlots, nPlotsPerFig) == 0
            fig = figure('Units', 'inches', 'Position', [0 0.5 20 10.0208], ...
                'Visible', 'off');  % 27 inch screen maximized window size
            tlo = tiledlayout(fig, nRow, nCol, 'TileSpacing', 'compact', 'Padding', 'compact');
        end

        % Plot vector field
        nexttile();
        lgd = PlotAttractorLandscape(sims(:, :, 1:nSimsShow), eq, lcMotif, 3, nSkipTR);
        fontsize(gca, basefont, "points");
        if nSess == 1
            title(sprintf('S%s, %dFP, %dLC', sublist(iSub), size(eq, 2), size(lc, 2)), ...
                'FontSize', basefont + 2);
        else
            title(sprintf('S%s, model %d, %dFP, %dLC', sublist(iSub), iSess, ...
                size(eq, 2), size(lc, 2)), 'FontSize', basefont + 2);
        end
        nPlots = nPlots + 1;

        % Close old figure
        if nPlots == nSubs * nSess || mod(nPlots, nPlotsPerFig) == 0
            lgd.Layout.Tile = 'east';
            lgd.Visible = 'on';
            title(tlo, sprintf('Vector field & motifs of %s', mdlName), ...
                'Interpreter', 'none', 'FontSize', basefont + 4, 'FontWeight', 'bold');
            disp('Saving figure ...');
            PrintAsSeen(fullfile(vfdir, ...
                sprintf('Traj_%s_%03d', mdlName, ceil(nPlots / nPlotsPerFig))), ...
                '-dpng', '-r200');
            close(fig);
        end      
    end
end


%% Save data

save(outFile, 'sublist', 'allEq', 'allLCMotif');

if ~vis
    return;
end


%% Show distribution of the number of FP and LC

nEq = cellfun(@(x) size(x, 2), allEq);
nLC = cellfun(@(x) size(x, 2), allLCMotif);
nEqLC = string(nEq(:)) + "FP, " + string(nLC(:)) + "LC";
nEqLC = categorical(nEqLC);

fig = figure('Units', 'inches', 'Position', [0 0.5 10 5], 'Visible', 'off');
histogram(nEqLC, 'Normalization', 'probability');
xlabel('Number of stable fixed points (FP) and stable limit cycles (LC)');
ylabel('Probability');
title(sprintf('Distribution of attractors in %s (n = %d)', mdlName, nSubs), ...
    'Interpreter', 'none');
PrintAsSeen(fullfile(figdir, sprintf('nFPLC_%s', mdlName)), '-dpng', '-r200');
close(fig);

end