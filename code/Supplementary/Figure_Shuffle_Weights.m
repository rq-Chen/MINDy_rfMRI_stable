function [tbl] = Figure_Shuffle_Weights(allMdl, method, p)
%FIGURE_SHUFFLE_WEIGHTS Get the distribution of dynamics in models with shuffled weights
%
%   Save down the number and type of attractors in the shuffled models.

[nSubs, nSess] = size(allMdl);
st = struct('Subject', [], 'Session', [], 'Method', [], 'Proportion', [], 'nFP', [], 'nLC', []);
st = repmat(st, nSubs * nSess, 1);
sMdl = GetShuffledModel(allMdl, method, p);
for i = 1:nSubs
    for j = 1:nSess
        idx = (i - 1) * nSess + j;
        st(idx).Subject = i;
        st(idx).Session = j;
        st(idx).Method = method;
        st(idx).Proportion = p;
        [eq, lc] = GetEqLC(sMdl{i, j});
        st(idx).nFP = size(eq, 2);
        st(idx).nLC = size(lc, 2);
    end
end
tbl = struct2table(st);

end