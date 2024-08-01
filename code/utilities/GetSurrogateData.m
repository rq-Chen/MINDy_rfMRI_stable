function [shuffled, matchSp, matchCov, matchMSpCov] = GetSurrogateData(dat, n)
%GETSURROGATEDATA Generate surrogate data
%
%   Input:
%       - |dat|: (nParcels, nTimepoint) data or cell array of it
%       - |n|: number of surrogate datasets to generate
%
%   Output:
%       - |shuffled|: temporally shuffled data of size (nParcels, nTimepoint, n)
%       - |matchSp|: simulated surrogate data of size (nParcels, nTimepoint,
%       n) that matches the spectrum of each parcel but with identity covariance
%       - |matchCov|: simulated surrogate data that matches the covariance
%       of the data but with white noise spectrum
%       - |matchSpCov|: simulated surrogate data that matches the AVERAGED
%       SPECTRUM of all parcels as well as the covariance of the data
%
%   Note that the function will linearly interpolate all the missing values
%   in |dat| before computation.
%
%   If |dat| is a cell array of data, it will be combined in time. The
%   results will be then splited in time back to cell arrays.

if iscell(dat)
    Ns = cellfun(@(x) size(x, 2), dat);
    dat = [dat{:}];
else
    Ns = [];
end
if nargin < 2; n = 1; end

% Handle input
dat = dat';
dat = fillmissing(dat, 'linear');
nParcels = size(dat, 2);
nTimepoint = size(dat, 1);

shuffled = nan(nParcels, nTimepoint, n);
for i = 1:n
    shuffled(:, :, i) = dat(randperm(nTimepoint), :)';
end
shuffled = SplitMat(shuffled, Ns);

if nargout > 1
    % FFT
    S = fft(dat);
    absS = abs(S);
    
    % White noise
    X = randn([size(dat) n]);
    
    % Match spectrum
    XS = fft(X);
    XSS = XS ./ abs(XS) .* absS;
    matchSp = ifft(XSS, "symmetric");
    matchSp = pagetranspose(matchSp);

    matchSp = SplitMat(matchSp, Ns);
end

% Match covariance
if nargout > 2
    C = cov(dat);
    [W, L] = eig(C);
    WsqL = W .* sqrt(diag(L)');
    matchCov = pagemtimes(WsqL, 'none', X, 'transpose');
    matchCov = SplitMat(matchCov, Ns);
end

% Match averaged spectrum and covariance
if nargout > 3
    XSS = XS .* mean(absS, 2);
    S = ifft(XSS, "symmetric");
    S = (S - mean(S)) ./ std(S);
    matchMSpCov = pagemtimes(WsqL, 'none', S, 'transpose');
    matchMSpCov = SplitMat(matchMSpCov, Ns);
end

end


%% Utilities
function dat = SplitMat(dat, Ns)
    if isempty(Ns)
        return
    end
    dat = mat2cell(dat, size(dat, 1), Ns(:), size(dat, 3));
    dat = reshape(dat, size(Ns));
end