function [allMdl, allDat, runIdx, sublist] = GetHCPRestModel(varargin)
%GETHCPRESTMODEL - Get deconvoluted data and MINDy models for HCP resting state fMRI
%
%   Usage: check GetHCPRestModel('help') for more information.
%
%   Default output location is "./data".

%% Parameters

p = inputParser();
p.FunctionName = 'GetHCPRestModel';
% Preprocessed data directory
addParameter(p, 'datDir', '../static/Singh2020PreProc/GSR3', @ischar);
% How the input was preprocessed, 'FIX' (Singh2020PreProc, with ICA-FIX) or 'NoFIX' (MINDyRestTaskPreProc)
addParameter(p, 'PreprocType', 'FIX', @(x) ismember(x, {'FIX', 'NoFIX'}));
% Subset of subject IDs to process (set as empty to use all available subjects)
addParameter(p, 'subListFile', '', @ischar);
% Subset of subject IDs to exclude (set as empty to exclude none)
addParameter(p, 'excludeListFile', '', @ischar);
% MINDy model type ('Simple', 'Masked' (MINDy_Rec), or 'HRF', only 'Simple' is supported for now)
addParameter(p, 'MINDyType', 'Simple', @(x) ismember(x, {'Simple', 'Masked', 'HRF'}));
% Number of parcels (100, 200, or 400)
addParameter(p, 'nParcels', 200, @(x) ismember(x, [100, 200, 400]));
% Which runs to use:
% - 'sess': two models for two sessions;
% - 'all': one model on all clean runs;
% - 'random': one model on |minRuns| (set as 2) random runs
addParameter(p, 'useRun', 'all', @(x) ismember(x, {'sess', 'all', 'random'}));
% Whether to output the BOLD data
addParameter(p, 'output_BOLD', false, @islogical);
% Whether to output the deconvoluted data
addParameter(p, 'output_Deconv', true, @islogical);
% Whether to output the MINDy models
addParameter(p, 'output_Mdl', true, @islogical);

% End of parameter definition

% Help
if nargin > 0 && strcmpi(varargin{1}, 'help')
    s = readlines([mfilename('fullpath'), '.m']);
    idx = find(strcmp(s, "% End of parameter definition"), 1, 'first');
    fprintf('%s\n', s{1:idx});
    return;
end

% Parse
parse(p, varargin{:});
datDir = p.Results.datDir;
PreprocType = p.Results.PreprocType;
subListFile = p.Results.subListFile;
excludeListFile = p.Results.excludeListFile;
MINDyType = p.Results.MINDyType;
nParcels = p.Results.nParcels;
useRun = p.Results.useRun;
output_BOLD = p.Results.output_BOLD;
output_Deconv = p.Results.output_Deconv;
output_Mdl = p.Results.output_Mdl;


% MINDy model
switch MINDyType
    case 'Simple'
        init_prj('MINDy_Base_v1.0');
    case 'Masked'
        init_prj('MINDy_RC');
    case 'HRF'
        init_prj('MINDy_HRF_v1.0');
end

% Training
nOrigRuns = 4;  % Number of runs in the original data (4 for HCP)
censorThres = 400;  % Exclude runs with this number of bad frames (total frames = 1200 for HCP)
if strcmpi(useRun, "sess")
    minRuns = nOrigRuns;  % Total number of (clean) runs needed (i.e., the threshold for excluding subjects)
elseif strcmpi(useRun, "random")
    minRuns = 2;
elseif strcmpi(useRun, "all")
    minRuns = 2;  % Models can be trained with different number of runs, but number of training batches are the same.
end
maskFile = fullfile('data', 'atlas', 'Wmask_RC.mat');  % Will be ignored unless |MINDyType| = 'Masked'
clearRecW = false;  % Clear parameter update history or not (only useful for 'Masked' model)

% Output
if output_BOLD
    preOut = fullfile('data', sprintf('HCP_Rest_%s_%s_BOLD%d_%s.mat', ...
        PreprocType, MINDyType, nParcels, useRun));  % BOLD data
else
    preOut = '';
end
if output_Deconv
    datOut = fullfile('data', sprintf('HCP_Rest_%s_%s_Deconv%d_%s.mat', ...
        PreprocType, MINDyType, nParcels, useRun));  % Deconvoluted data
else
    datOut = '';
end
if output_Mdl
    mdlOut = fullfile('data', sprintf('HCP_Rest_%s_%s_Mdl%d_%s.mat', ...
        PreprocType, MINDyType, nParcels, useRun));  % Models
else
    mdlOut = '';
end

% Utilify function
dropSub = @(xx)(xx(1:nParcels,:));  % Drop last 19 rows of subcortical regions


%% Get subject list

if strcmpi(PreprocType, "FIX")
    suffix = sprintf("Y0%d.mat", nParcels / 100);
elseif strcmpi(PreprocType, "NoFIX")
    suffix = sprintf("Y0%d_All.mat", nParcels / 100);
end

if ~isempty(subListFile)
    sublist = "sub" + readlines(subListFile) + suffix;  % string array
else
    tmp = dir(fullfile(datDir, "sub*" + suffix));
    sublist = string({tmp.name})';
end

if ~isempty(excludeListFile)
    excludeList = "sub" + readlines(excludeListFile) + suffix;
    sublist = setdiff(sublist, excludeList);
end

nSubs = numel(sublist);

InputFiles = fullfile(datDir, sublist);


%% Load and handle data

PreDat = cell(nSubs, nOrigRuns);
Ncensored = nan(nSubs, nOrigRuns);
TMask = cell(nSubs, nOrigRuns);
for i = 1:nSubs
    load(InputFiles(i), 'X');
    disp("Loading " + sublist(i));
    if strcmpi(PreprocType, "NoFIX")
        X.Dat = X.Dat{1};
        X.QC = X.QC{1};
        if isempty(X.Dat) || isempty(X.QC)
            continue;
        end
    end
    tmpIdx = find(cellfun(@(x) ~isempty(x), X.Dat));  % Use find, since length(X) = 4 but length(X.QC.run) <= 4
    Ncensored(i, tmpIdx) = arrayfun(@(x) sum(x.tmask), X.QC.run(tmpIdx));
    TMask(i, tmpIdx) = arrayfun(@(x) x.tmask, X.QC.run(tmpIdx), "UniformOutput", false);
    PreDat(i, tmpIdx) = cellfun(dropSub, X.Dat(tmpIdx), 'UniformOutput', false);
end

% Remove bad runs altogether (bad frames will be handled in MINDy fitting)
is_good_run = Ncensored < censorThres;  % Note: this will be false for Ncensored == nan, so OK
PreDat(~is_good_run) = {[]};
% Remove subjects without enough runs
tmpIdx = sum(is_good_run, 2) >= minRuns;
Ncensored = Ncensored(tmpIdx, :);
TMask = TMask(tmpIdx, :);
PreDat = PreDat(tmpIdx, :);
InputFiles = InputFiles(tmpIdx);
sublist = sublist(tmpIdx);
is_good_run = is_good_run(tmpIdx, :);
nSubs = length(sublist);

disp("Number of subjects: " + string(nSubs));
if nSubs > 250
    warning('More than 250 subjects. Data file will be large (>1GB for 100 parcel model).')
end

% Save BOLD data
if ~isempty(preOut)
    save(preOut, 'PreDat', 'Ncensored', 'TMask', 'InputFiles', 'sublist', '-v7.3');
end


%% Training

% Wmask
if nParcels == 100 && ~isempty(maskFile)
    load(maskFile, 'Wmask');
else
    Wmask = 1;
end

% Run index
if strcmp(useRun, "random")
    runIdx = cell(nSubs, 1);
    for i = 1:nSubs
        tmp = find(is_good_run(i, :));
        runIdx{i} = tmp(randperm(length(tmp), minRuns));
    end
elseif strcmp(useRun, "sess")
    runIdx = repmat({1:2, 3:4}, nSubs, 1);
elseif strcmp(useRun, "all")
    runIdx = cell(nSubs, 1);
    for i = 1:nSubs
        runIdx{i} = find(is_good_run(i, :));
    end
end

% Training
[allMdl, allDat] = GetMINDy(PreDat, runIdx, Wmask, clearRecW, MINDyType);

% Output
if ~isempty(datOut)
    save(datOut, "allDat", "runIdx", "InputFiles", "sublist", '-v7.3');
end
if ~isempty(mdlOut)
    save(mdlOut, "allMdl", "runIdx", "InputFiles", "sublist", '-v7.3');
end

end