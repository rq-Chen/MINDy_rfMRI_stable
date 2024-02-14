function Dat = MINDy_SimpleTaskPreProcInterp(Dat,TR,zscoring,nEdgeTRs,varargin)
%MINDY_SIMPLETASKPREPROCINTERP Preprocessing tfMRI data in a way similar to MINDy_Simple() does
%
%   Inputs:
%       Dat: (nParcels, nTRs) matrix or cell array of them.
%       TR: TR of the data in seconds, default is set in ChosenPARSTR(), which is 0.72.
%       zscoring: whether to zscore the data, default is False.
%       nEdgeTRs: number of TRs to discard at both ends of the data, default is half of the HRF length.
%
%   Output:
%       Dat: preprocessed data, same size as the input.
%
%   Pipeline (similar to MINDy_Simple() except that NaNs are not discarded):
%       1. Interpolate NaNs (NO extrapolation). The indices of remaining NaNs at both ends are recorded.
%       2. zscoring (optional).
%       3. Call MINDy_RestingStatePreProcInterp() WITHOUT extreme value filtering. Will return de-conved
%          data of the same length as the non-NaN part of the original data.
%       4. Set the first and last |nEdgeTRs| TRs of the returned data to NaNs.
%       5. Insert NaNs back at the indices recorded in step 1.

% Some parameters
ChosenPARSTR;
Pre.FiltAmp = inf;  % No extreme value filtering
HRF_length = 30;  % In TRs

% Handle inputs
if nargin < 4 || isempty(nEdgeTRs)
    nEdgeTRs = round(HRF_length / 2);
end
if nargin < 3 || isempty(zscoring)
    zscoring = false;
end
if nargin < 2 || isempty(TR)
    TR = Pre.TR;
else
    Pre.TR = TR;
end
if ~iscell(Dat)
    not_cell_input = true;
    Dat = {Dat};
else
    not_cell_input = false;
end

% Preprocessing
for i = 1:numel(Dat)

    dat = Dat{i};  % (nParcels, nTRs)

    % Interpolate NaNs
    for j = 1:size(dat, 1)
        nanIdx = isnan(dat(j, :));
        dat(j, nanIdx) = interp1(find(~nanIdx), dat(j, ~nanIdx), find(nanIdx));
    end

    % Get the index of remaining frames that contain NaNs
    nanIdx = any(isnan(dat));
    dat(:, nanIdx) = nan;

    % zscoring
    if zscoring
        mu = mean(dat, 2, "omitnan");
        sigma = std(dat, 0, 2, "omitnan");
        dat = (dat - mu) ./ sigma;
    end

    % Preprocessing
    newDat = MINDy_RestingPreProcInterp({dat(:, ~nanIdx)}, Pre.FiltAmp, Pre.ConvLevel, Pre.DownSamp, TR, HRF_length);
    newDat = newDat{1};

    % Discard edge frames
    newDat(:, 1:nEdgeTRs) = nan;
    newDat(:, end - nEdgeTRs + 1:end) = nan;

    % Insert NaNs back
    dat(:, ~nanIdx) = newDat;
    Dat{i} = dat;

end

if not_cell_input
    Dat = Dat{1};
end

end