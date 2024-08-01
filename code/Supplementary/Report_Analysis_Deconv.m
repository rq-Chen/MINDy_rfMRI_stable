function [] = Report_Analysis_Deconv(mdlName, mdlFile, datFile, vis_dir, varargin)
%REPORT_ANALYSIS_DECONV Some analysis directly on data
%
%   Report_Analysis_Deconv(mdlName, mdlFile, datFile, vis_dir, varargin)
%
% Positional input:
%   mdlName: str, the name of the model
%   mdlFile: str, the path to the model file
%   datFile: str, the path to the data file
%   vis_dir: str, the directory to save the figures


%% Parameters

if nargin < 1 || isempty(mdlName)
    mdlName = 'HCP_Rest_FIX_Simple_Mdl200_sess';
end
if nargin < 2 || isempty(mdlFile)
    mdlFile = fullfile('data', [mdlName '.mat']);
end
if nargin < 3 || isempty(datFile)
    datFile = fullfile('data', [replace(mdlName, 'Mdl', 'Deconv') '.mat']);
end
if nargin < 4 || isempty(vis_dir)
    vis_dir = fullfile('figures', mdlName, 'FC');
end
if ~exist(vis_dir, "dir")
    mkdir(vis_dir);
end

% Number of participants to process
nSubs = 30;
nSess = 2;


%% Load models & data

load(mdlFile, 'allMdl', 'runIdx');
load(datFile, 'allDat');

% Calculate RMSE as noise std
if ~isfield(allMdl{1}, 'RMSE')
    for i = 1:size(runIdx, 1)
        for j = 1:size(runIdx, 2)
            xt = cellfun(@(x) x(:, 1:end-2), allDat(i, runIdx{i, j}), 'UniformOutput', false);
            xt = [xt{:}];
%             dxt = cellfun(@(x) x(:, 2:end-1) - x(:, 1:end-2), allDat(i, runIdx{i, j}), 'UniformOutput', false);
            dxt = cellfun(@(x) (x(:, 3:end) - x(:, 1:end-2)) / 2, allDat(i, runIdx{i, j}), 'UniformOutput', false);
            dxt = [dxt{:}];
            pred = MINDyInt_00(allMdl{i, j}, xt, 1, 1, 0, 1);
            pred = squeeze(pred(:, 2, :) - pred(:, 1, :));
            allMdl{i, j}.RMSE = sqrt(mean((dxt - pred) .^ 2, 2));
        end
    end
    save(mdlFile, 'allMdl', '-append');
end


%% Simulate models with noise

allDat = allDat(1:nSubs, :);
runIdx = runIdx(1:nSubs, :);
allMdl = allMdl(1:nSubs, :);

nInits = 2;  % Number of runs to simulate per model
% option = optimset('Display', 'iter', 'PlotFcns', {@optimplotx, @optimplotfval}, ...
%     'TolX', 1e-3, 'MaxFunEvals', 50);

[noise_scale, simDat, simFC, datFC] = GetOptimalNoise(allMdl, ...
    allDat, runIdx, nInits, 0.1, 1, true);
simDat = reshape(permute(simDat, [1 3 2]), nSubs, []);
fprintf('Noise scaling factor = %.3f\n', noise_scale);


%% FC correlation between data and simulated data

trilv = @(x) squareform(tril(x, -1))';

FCcorr = cellfun(@(x, y) corr(trilv(x), trilv(y)), simFC, datFC);

% Show results
histogram(FCcorr(:), 'Normalization', 'probability');
xlabel('Correlation between data and simulated FCs');
ylabel('Proportion of subjects-sessions');
title('FC Correlation between Data and Simulated Data');
subtitle(sprintf('Subset of %d subjects', nSubs));
saveas(gcf, fullfile(vis_dir, 'FC_corr_histogram.png'));
close(gcf);


%% Dynamic FC state analysis per subject-session

allK = 1:5;  % List of the number of clusters
winlen = 52;
lag = 52;

C1 = dFC_state_analysis(allDat, allK, winlen, lag);
C2 = dFC_state_analysis(simDat, allK, winlen, lag);
PlotDFCStates(C1, C2);

PrintAsSeen(fullfile(vis_dir, sprintf('dFCstates_winlen-%d_lag-%d', winlen, lag)), '-dpng', '-r300');
PrintAsSeen(fullfile(vis_dir, sprintf('dFCstates_winlen-%d_lag-%d', winlen, lag)), '-dpdf', '-vector');
close(gcf);


end