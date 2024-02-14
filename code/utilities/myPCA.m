function [score, coeff, explained, mu] = myPCA(sims, K, coeff, m)
% MYPCA Compute PCA for simulations
%   Input:
%       |sims|: (nParcels, T, nSims) or cell array of it (Note: nSims can be 1)
%       |K|: number of components required, default is all
%       |coeff|: optional, use |coeff| as PC coefficients instead
%       |m|: optional, use |m| as the mean instead (only used with |coeff|)
%   Output:
%       |score|: (T, K) if nSims == 1, or (T, nSims, K), or cell array of it
%       |coeff|: (nParcels, K)
%       |explained|: (nParcels, 1) or (length(coeff), 1) if |coeff| is provided
%       |mu|: (1, nParcels);
    if nargin == 1 || (nargin == 2 && isempty(K))
        if iscell(sims)
            K = size(sims{1}, 1);
        else
            K = size(sims, 1);
        end
    elseif isempty(K)
        K = size(coeff, 2);
    end

    if iscell(sims)
        nSimsvec = cellfun(@(x) size(x, 3), sims);
        tmp = permute(cat(3, sims{:}), [2 3 1]);
    else
        tmp = permute(sims, [2 3 1]);
    end
    [T, SS, nParcels] = size(tmp);
    tmp = reshape(tmp, [], nParcels);
    if nargin >= 3
        if nargin >= 4
            mu = m;
            if iscolumn(mu)
                mu = mu';
            end
        else
            mu = mean(tmp);
        end
        coeff = coeff(:, 1:K);
        score = (tmp - mu) * coeff;
        datCov = cov(tmp);
        explained = diag(coeff' * datCov * coeff) / trace(datCov) * 100;
    else
        warning('off', 'stats:pca:ColRankDefX');
        [coeff, score, ~, ~, explained, mu] = pca(tmp, ...
            'NumComponents', K);
        warning('on', 'stats:pca:ColRankDefX');
    end
    score = reshape(score, T, SS, K);
    if SS == 1
        score = squeeze(score);
    end
    if iscell(sims)
        score = mat2cell(score, T, nSimsvec, K);
    end
end