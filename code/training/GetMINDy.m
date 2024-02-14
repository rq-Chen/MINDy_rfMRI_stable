function [allMdl, allDat] = GetMINDy(PreDat, runIdx, Wmask, clearRecW, MINDyType)
%GETMINDY Train MINDy models from parcellated data
%
%   Update 2023/08/23: now |runIdx| can be a (nSubs, N) cell array, so you can
%   use different runs and even different number of runs for different subjects.
%   The new default for |runIdx| is to use all nonempty cells in |PreDat|.
%
%   This function trains required type |MINDyType| of MINDy models using the
%   data in |PreDat|. It trains |N| models for each subject where |N| is the
%   number of cells in |runIdx|. |runIdx{i}| specifies which runs to use to
%   train the |i|-th model (default: use all). The function outputs all trained
%   models as a cell array and the processed data (in neural space) as another
%   cell array.
%
%   Inputs:
%       - |PreDat|: (nSubs, nRuns) cell array of (nParcels, nTimepoints) matrices
%       - |runIdx|: (1, N) cell array of index vectors. Example: {[1 2], [3 4]}
%           for two models per subject, each trained on one session.
%       - |Wmask|: sparsity mask to the connectivity matrix. Default is no mask.
%           Only used when |MINDyType| is 'Masked'.
%       - |clearRecW|: whether to clear the parameter update history to save space.
%       - |MINDyType|: 'Simple', 'Masked' or 'HRF' (default: 'Simple').
%
%   Outputs:
%       - |allMdl|: (nSubs, N) cell array of MINDy model structures
%       - |allDat|: (nSubs, nRuns) cell array of processed data. Same order as
%           |PreDat|.

    if nargin < 5
        MINDyType = 'Simple';
    end
    if nargin < 4
        clearRecW = true;
    end
    if nargin < 3 || isempty(Wmask)
        Wmask = 1;
    end
    if nargin < 2 || isempty(runIdx)
        runIdx = cell(size(PreDat, 1), 1);
        for i = 1:size(PreDat, 1)
            runIdx{i} = find(~cellfun(@isempty, PreDat(i, :)));
        end
    end
    if size(runIdx, 1) == 1
        runIdx = repmat(runIdx, size(PreDat, 1), 1);
    end

    nSub = size(runIdx, 1);
    nSess = size(runIdx, 2);
    % nParcels = size(PreDat{1}, 1);

    allMdl = cell(nSub, nSess);
    if nargout > 1
        allDat = cell(size(PreDat));
    end
    
    for i = 1:nSub
        for j = 1:nSess
            if strcmpi(MINDyType, 'Masked')
                [allMdl{i, j}, ~, ~, tmpDat] = ...
                    MINDy_Simple_Rec(PreDat(i, runIdx{i, j}), .72, 'n', Wmask, nBatch);
            elseif strcmpi(MINDyType, 'Simple')
                [allMdl{i, j}, ~, ~, tmpDat] = ...
                    MINDy_Simple(PreDat(i, runIdx{i, j}), .72, 'y');
            elseif strcmpi(MINDyType, 'HRF')
                [tmpMdl, ~, ~, ~, tmpDat] = ...
                    MINDy_HRFbold_Simple_Din(PreDat(i, runIdx{i, j}), .72, 'y');
                % Remove data and HRFout (data will be saved separately; HRFout is not used but HUGE)
                allMdl{i, j} = rmfield(tmpMdl, {'deconvX', 'HRFout'});
            end
            % Save processed data
            if nargout > 1
                allDat(i, runIdx{i, j}) = tmpDat;  % Note: tmpDat is a cell array
            end
            % Clear up unused field to save memory and parpool overhead
            if clearRecW && isfield(allMdl{1}, 'RecW')
                allMdl{i, j} = rmfield(allMdl{i, j}, {'RecW', 'RecA', 'RecD'});
            end
        end
    end
end


