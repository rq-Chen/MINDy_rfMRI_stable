function [types, subjID, iSess] = get_example_model_idx(plotType)
    if nargin == 0 || strcmp(plotType, 'popular')
        types = {'2FP 0LC', '0FP 1LC', '1FP 0LC', '4FP 0LC'};
        subjID = ["100206" "108323" "100307" "110411"];
        iSess = [1 2 2 2];
    elseif strcmp(plotType, 'rare')
        types = {'0FP 2LC', '2FP 2LC', '2FP 1LC', '6FP 0LC'};
        subjID = ["618952" "105923" "734045" "157437"];
        iSess = [1 1 2 2];
    end
end