%% Load data
function [allDat, allMdl, runIdx] = LoadData(datFile, mdlFile, clearRecW)

    if nargin < 3
        clearRecW = true;
    end

    if ~isempty(datFile)
        load(datFile, 'allDat');
    else
        allDat = [];
    end

    if ~isempty(mdlFile)
        warning('off', 'MATLAB:load:variableNotFound');
        tmp = load(mdlFile, 'Val', 'runIdx', 'allMdl');
        warning('on', 'MATLAB:load:variableNotFound');
        if isfield(tmp, 'runIdx')
            runIdx = tmp.runIdx;
        else
            runIdx = {[1 2], [3 4]};
        end
        if isfield(tmp, 'Val')  
            if isstruct(tmp.Val)
                allMdl = tmp.Val.Val;
            else
                allMdl = tmp.Val;
            end
        else
            allMdl = tmp.allMdl;
        end
        
        % Clear up unused field in |allMdl| to save memory and parpool overhead
        if clearRecW && isfield(allMdl{1}, 'RecW')
            allMdl = cellfun(@(x) rmfield(x, {'RecW', 'RecA', 'RecD'}), allMdl, ...
                'UniformOutput', false);
        end
    end
end


