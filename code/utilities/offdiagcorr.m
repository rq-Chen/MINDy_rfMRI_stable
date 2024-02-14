function r = offdiagcorr(A, B)
assert(isequal(size(A), size(B)), 'A and B must have the same size!');
for i = 1:size(A, 1)
    A(i, i) = nan;
    B(i, i) = nan;
end
r = corr(A(:), B(:), 'rows', 'complete');
end