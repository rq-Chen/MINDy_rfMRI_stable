function im = SurfacePlotIC(lmst, nFrames, select_by, start_end, t, mmpL, mmpR, sL, sR)
%%SURFACEPLOTIC Surface plot (gif) for invariant cycles
%
% WARNING: FLICKERING when generating plots.
%
% Inputs:
%   - lmst: Limit set cell or structure
%   - nFrames: number of frames
%   - select_by: "time" or "distance"
%   - start_end: "start" or "end", used to show two directions from a saddle
%   - t: title prefix

if iscell(lmst)
    dat = lmst{1};
    tmp = 0:(nFrames - 1); 
    Idx = round(size(dat, 2) / nFrames * tmp) + 1;
elseif isstruct(lmst)
    looplen = cellfun(@length, lmst.C);
    [icLen, idx] = max(looplen);
    icV = lmst.C{idx};

    % The fixed points on the cycle
    icE = nan(icLen, 2);
    icE(:, 1) = icV;
    icE(1:end - 1, 2) = icV(2:end);
    icE(end, 2) = icV(1);
    % Flip the start and end so saddle is always the first
    if lmst.VUSDim(icE(1, 1)) == 0  % Att-Saddle-Att...
        icE(1:2:end, :) = icE(1:2:end, [2 1]);
    else  % Saddle-Att-Saddle...
        icE(2:2:end, :) = icE(2:2:end, [2 1]);
    end

    if isempty(start_end) || strcmpi(start_end, 'start')
        dat = lmst.E{icE(1, 1), icE(1, 2)};
    else
        if icE(end, 1) == icE(1, 1)
            dat = lmst.E{icE(end, 1), icE(end, 2)};
        else
            dat = lmst.E{icE(2, 1), icE(2, 2)};
        end
    end
    
    % (1, nTimepoint) traveled distance
    dist = cumsum(vecnorm([dat(:, 2:end) - dat(:, 1:end - 1) zeros(size(dat, 1), 1)]));
    
    % Select frames by percentile of the distance or time
    if strcmpi(select_by, 'distance')
        q = dist(end) * linspace(0, 1, nFrames);
        Idx = nan(size(q));
        Idx(1) = 1; Idx(end) = length(q);
        tmpI = 2;
        for i = 2:length(dist)
            if (dist(i) > q(tmpI))
                Idx(tmpI) = i;
                tmpI = tmpI + 1;
                if tmpI > length(q) - 1
                    break
                end
            end
        end
    elseif strcmpi(select_by, 'time')
        % Discard the (1/nFrames) distance at both ends as they are too slow
        q = dist(end) * [1 / nFrames 1 - 1 / nFrames];
        [~, qIdx] = min(abs(dist' - q));
        Idx = [1 round(linspace(qIdx(1), qIdx(2), nFrames - 2)) length(dist)];
    end
end

dat = dat(:, Idx);

% Plot
im = cell(nFrames, 1);
for i = 1:nFrames
    f = figure(); f.WindowState = 'maximized';
    PlotYeoSurface(dat(:, i), mmpL, mmpR, sL, sR);
    clim([-1 1]); colorbar;
    title(sprintf('%s, TR = %d', t, Idx(i)))
    drawnow;
    im{i} = frame2im(getframe(f));
    close(f);
end

end