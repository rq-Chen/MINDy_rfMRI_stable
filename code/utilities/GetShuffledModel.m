function [sMdl] = GetShuffledModel(allMdl, method, p)
%GETSHUFFLEDMODEL Get model with shuffled weights
%
%    [sMdl] = GetShuffledModel(allMdl, method, p);
%
% Inputs:
%    allMdl: model structure or cell array of model structures
%    method: "within_model", "within_person" or "across"
%    p: proportion of weights to shuffle
%
% Outputs:
%    sMdl: model structure with shuffled weights, or cell array of models
%
% For "within_model" method, the weights of each model are shuffled across the
% entries of the weight matrix. For "within_person" method, the weights will
% be shuffled across the corresponding entries of the weight matrices of the
% models from the same person (i.e., in the same row). For "across" method,
% the weights were shuffled across the corresponding entries of the weight
% matrices of all models.

if isstruct(allMdl)
    allMdl = {allMdl};
end
if nargin < 2
    method = 'within_model';
end
if nargin < 3
    p = 1;
end
paramIdx = [2 5 6];  % Indices of fields in mdl.Param to shuffle
shuffleMask = cellfun(@(x) rand(size(x)) < p, allMdl{1}.Param(paramIdx), 'UniformOutput', false);
shuffleIdx = cellfun(@(x) find(x), shuffleMask, 'UniformOutput', false);
[clIdx, rIdx] = cellfun(@(x) find(x), shuffleMask, 'UniformOutput', false);

sMdl = cell(size(allMdl));
if strcmp(method, 'within_model')
    for i = 1:numel(allMdl)
        tmp = struct();
        tmp.Param = allMdl{i}.Param;
        param = allMdl{i}.Param(paramIdx);
        for j = 1:numel(param)
            perm = randperm(numel(shuffleIdx{j}));
            param{j}(shuffleIdx{j}) = param{j}(shuffleIdx{j}(perm));
        end
        tmp.Param(paramIdx) = param;
        sMdl{i} = tmp;
    end
elseif strcmp(method, 'within_person')
    for i = 1:size(allMdl, 1)
        % Gather parameters from all models of the same person
        params = cell(1, numel(paramIdx));
        for j = 1:numel(paramIdx)
            params{j} = cellfun(@(x) x.Param{paramIdx(j)}, allMdl(i, :), 'UniformOutput', false);
            params{j} = cat(3, params{j}{:});
            for l = 1:numel(clIdx{j})
                perm = randperm(size(params{j}, 3));
                params{j}(clIdx{j}(l), rIdx{j}(l), :) = params{j}(clIdx{j}(l), rIdx{j}(l), perm);
            end
        end
        % Assign shuffled parameters to models
        for j = 1:size(allMdl, 2)
            S = struct();
            S.Param = allMdl{i, j}.Param;
            for k = 1:numel(paramIdx)
                S.Param{paramIdx(k)} = params{k}(:, :, j);
            end
            sMdl{i, j} = S;
        end
    end
elseif strcmp(method, 'across')
    % Gather parameters from all models
    params = cell(1, numel(paramIdx));
    for j = 1:numel(paramIdx)
        params{j} = cellfun(@(x) x.Param{paramIdx(j)}, allMdl, 'UniformOutput', false);
        params{j} = cat(3, params{j}{:});
        for l = 1:numel(clIdx{j})
            perm = randperm(size(params{j}, 3));
            params{j}(clIdx{j}(l), rIdx{j}(l), :) = params{j}(clIdx{j}(l), rIdx{j}(l), perm);
        end
    end
    % Assign shuffled parameters to models
    for i = 1:numel(allMdl)
        S = struct();
        S.Param = allMdl{i}.Param;
        for j = 1:numel(paramIdx)
            S.Param{paramIdx(j)} = params{j}(:, :, i);
        end
        sMdl{i} = S;
    end
else
    error('Unknown method: %s', method);
end

end