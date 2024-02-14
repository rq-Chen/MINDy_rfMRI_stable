function [] = GetSubjMotifs(mdlName, varargin)
%GETSUBJMOTIFS Save subject ID and individualized motifs in a file

%% Parameters

if nargin == 0 || isempty(mdlName)
    mdlName = 'HCP_Rest_FIX_Simple_Mdl200_sess';
end

p = inputParser;
addParameter(p, 'mtfFile', fullfile('data', ['motifs_' mdlName '.mat']), @ischar);
addParameter(p, 'outFile', fullfile('data', ['subjMotifs_' mdlName '.mat']), @ischar);
addParameter(p, 'selection', 'max_sim', @(x) ismember(x, {'max_sim', 'min_sim', 'mean'}))

parse(p, varargin{:});
mtfFile = p.Results.mtfFile;
outFile = p.Results.outFile;
selection = p.Results.selection;


%% Calculation

% Get motifs and do clustering
[V, ~, clIdx, C] = Report_Motifs(mdlName, 'mtfFile', mtfFile, 'vis', false, ...
    'normalize', false, 'exclude_origin', true, 'mtfType', 'all');
[K, nParcels] = size(C);
[nSubs, nSess] = size(V);

% Load subject IDs
load(mtfFile, 'sublist');  % (nSubs, 1) string
assert(nSubs == length(sublist));

% Use the motif in the cluster but least similar to the centroid for each subject
sessMotifs = nan(nParcels, K, nSubs, nSess);
subjMotifs = nan(nParcels, K, nSubs);
for i = 1:nSubs
    for k = 1:K
        if strcmp(selection, 'mean')
            for j = 1:nSess
                sessMotifs(:, k, i, j) = mean(V{i, j}(:, clIdx{i, j} == k), 2);
            end
            tmpV = [V{i, :}];
            tmpIdx = [clIdx{i, :}];
            subjMotifs(:, k, i) = mean(tmpV(:, tmpIdx == k), 2);
        else
            optDist = nan;
            optIdx = [nan nan];  % Session, motif
            for j = 1:nSess
                idx = find(clIdx{i, j} == k);
                if isempty(idx)
                    continue;
                end
                distMat = pdist2(V{i, j}(:, idx)', C(k, :), "cosine");
                if strcmp(selection, 'max_sim')
                    [mDist, mIdx] = max(distMat);
                    mIdx = idx(mIdx);
                    if isnan(optDist) || mDist > optDist
                        optDist = mDist;
                        optIdx = [j mIdx];
                    end
                elseif strcmp(selection, 'min_sim')
                    [mDist, mIdx] = min(distMat);
                    mIdx = idx(mIdx);
                    if isnan(optDist) || mDist < optDist
                        optDist = mDist;
                        optIdx = [j mIdx];
                    end
                end
                sessMotifs(:, k, i, j) = V{i, j}(:, mIdx);
            end
            if ~isnan(optIdx(1))
                subjMotifs(:, k, i) = V{i, optIdx(1)}(:, optIdx(2));
            end
        end
    end
end

save(outFile, 'subjMotifs', 'sessMotifs', 'sublist', 'C');

end