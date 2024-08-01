%% Initial conditions sampled from data
function allInits = GetInits(allDat, nSims, runIdx, sim_per, reuse_init)

    % Inputs
    if nargin < 5
        reuse_init = true;
    end
    if nargin < 4
        sim_per = 'subject';
    end
    if nargin < 3
        runIdx = {[1 2], [3 4]};
    end
    if nargin < 2
        nSims = 120;
    end

    nSub = size(allDat, 1);
    nSess = size(runIdx, 2);
    if size(runIdx, 1) == 1
        runIdx = repmat(runIdx, nSub, 1);
    end

    allInits = cell(nSub, nSess);

    switch sim_per
        case 'population'
            currDat = cat(2, allDat{:});
            currL = size(currDat, 2);
            for i = 1:nSub
                for j = 1:nSess
                    if reuse_init && i * j > 1
                        allInits{i, j} = allInits{1, 1};
                    else
                        allInits{i, j} = currDat(:, randperm(currL, min(currL, nSims)));
                    end
                end
            end
        case 'subject'
            for i = 1:nSub
                currDat = cat(2, allDat{i, :});
                currL = size(currDat, 2);
                for j = 1:nSess
                    if reuse_init && j > 1
                        allInits{i, j} = allInits{i, 1};
                    else
                        allInits{i, j} = currDat(:, randperm(currL, min(currL, nSims)));
                    end
                end
            end
        case 'model'
            for i = 1:nSub
                for j = 1:nSess
                    currDat = cat(2, allDat{i, runIdx{i, j}});
                    currL = size(currDat, 2);
                    allInits{i, j} = currDat(:, randperm(currL, min(currL, nSims)));
                end
            end
    end
end


