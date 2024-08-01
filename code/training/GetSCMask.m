%GETSCMASK Generate a mask for W matrix based on structural connectivity
%
% [SC, Wmask, tSC, pSC] = GetSCMask(nParcels, MICA_dir, pt, WmaskFile)
%
% Inputs:
% - nParcels: number of parcels in the parcellation (default: 200)
% - MICA_dir: directory containing the MICA-MICs dataset (default: '../micapipe/')
% - pt: a percentile for thresholding the SC matrix (e.g. 90) or 'tstat'
% which thresholds the matrix based on Bonferroni-corrected t-statistic of
% log-transformed SC.
% - WmaskFile: file to save the mask (default: 'data/atlas/MICA{pt}.mat')
%
% Outputs:
% - SC: mean structural connectivity matrix (# of streamlines)
% - Wmask: mask for the W matrix
% - tSC: t-statistics of the log-transformed SC matrix
% - pSC: p-values of the t-statistics of the log-transformed SC matrix (uncorrected,
%        one-tailed)
%
% The mask is generated based on the mean SC matrix from the MICA-MICs dataset
% at https://osf.io/x7qr2 using the micapipe toolbox. The lookup table for the
% parcellation is available at "micapipe/parcellations/lut/" in micapipe's
% repository.
%
% The ordering of the rows and columns is described in the associated paper
% https://www.nature.com/articles/s41597-022-01682-y#Sec27, which is:
% - 7 left hemisphere (LH) subcortical regions
% - 7 right hemisphere (RH) subcortical regions
% - 1 medial wall
% - LH cortical regions
% - 1 medial wall
% - RH cortical regions

function [SC, Wmask, tSC, pSC] = GetSCMask(nParcels, MICA_dir, pt, WmaskFile)

if nargin < 1 || isempty(nParcels)
    nParcels = 200;
end
if nargin < 2 || isempty(MICA_dir)
    MICA_dir = '../micapipe/';
end
if nargin < 3 || isempty(pt)
    pt = 'tstat';
end
if nargin < 4 || isempty(WmaskFile)
    if isnumeric(pt)
        WmaskFile = fullfile('data', 'atlas', sprintf('MICA%02d.mat', pt));
    else
        WmaskFile = fullfile('data', 'atlas', 'MICAtstat.mat');
    end
end

% Load individual SC matrices
subjects = getFolders(MICA_dir);
nSubs = length(subjects);
allSC = cell(nSubs, 1);
for i = 1:nSubs
    subj = subjects(i);
    sessions = getFolders(fullfile(MICA_dir, subj));
    SCs = cell(1, length(sessions));
    for j = 1:length(sessions)
        sess = sessions(j);
        filename = subj + "_" + sess + "_space-dwinative_atlas-schaefer" + nParcels + "_desc-SC.txt";
        SCs{j} = readmatrix(fullfile(MICA_dir, subj, sess, "dwi", filename));
    end
    SCs = cat(3, SCs{:});
    SCs = mean(SCs, 3);
    allSC{i} = SCs;
end
allSC = cat(3, allSC{:});
SC = mean(allSC, 3);
tSC = mean(log(allSC), 3) ./ (std(log(allSC), [], 3) + eps) * sqrt(nSubs);

% Remove subcortical regions and medial wall
tmpIdx = nParcels / 2 + 15;
SC = SC([16:tmpIdx tmpIdx+2:end], [16:tmpIdx tmpIdx+2:end]);
tSC = tSC([16:tmpIdx tmpIdx+2:end], [16:tmpIdx tmpIdx+2:end]);

% Recover the lower triangle
SC = triu(SC) + tril(SC', -1);
tSC = triu(tSC) + tril(tSC', -1);

% Get p-value
pSC = tcdf(tSC, nSubs - 1, "upper");

% Threshold the SC matrix to obtain mask
if isnumeric(pt)
    thres = prctile(SC(:), pt);
    Wmask = (SC >= thres);
elseif ischar(pt) && strcmp(pt, 'tstat')
    thres = 0.05 / sum(triu(true(size(SC))), "all");
    Wmask = (pSC < thres);
end

% Save mask
save(WmaskFile, "Wmask");

end


%% Utility functions

function [folders] = getFolders(path)
    folders = dir(path);
    folders = folders([folders.isdir]);
    folders = string({folders.name});
    folders = folders(~ismember(folders, ["." ".."]));
end