function motifs = LC2Motif(LC, motif_type, varargin)
%LC2MOTIF Converts limit cycles to their motifs
%
%   Update: move the center from the first output to the last
%     and omit it by default.
%
%   Usage:
%
%     motifs = LC2Motif(LC, motif_type);
%     motifs_w_center = LC2Motif(LC, motif_type, "w_center");
%
%   Input:
%
%     LC: limit cycles. (1, nLC) cell array of (nParcels, 
%       tPeriod) matrices.
%     motif_type: string or character array. Method to get motifs.
%       - "unscaled_centered_PC"
%       - "unscaled_centered_combined_PC"
%       - "PC_max_projection"
%       - "Point_max_projection"
%       - "PC_max_cosine"
%       - "Point_max_cosine"
%       - "Slowest" (default)
%       It's easier to understand by reading the codes.
%
%   Output:
%
%     motifs: motifs. (1, nLC) cell array of (nParcels, 4)
%       matrices. The columns are:
%         1. "positive semi-major extreme"
%         2. "positive semi-minor extreme" 
%         3. "negative semi-major extreme"
%         4. "negative semi-minor extreme"
%
%     motifs_w_center: similar to |motifs| but with an extra column:
%         5. "center" of the limit cycle;

% Handle inputs
assert(iscell(LC), "LC should be a cell array of matrices.")
if nargin < 2 || isempty(motif_type)
    motif_type = "Slowest";
end
motif_type = string(motif_type);
if strcmpi(motif_type, "unscaled_centered_PC") && numel(LC) > 1
    warning("unscaled_centered_PC is not suitable for LCs that are not centered at the origin.");
end
if strcmpi(motif_type, "unscaled_centered_combined_PC")
    LC = {[LC{:}]};
    motif_type = "unscaled_centered_PC";
end
if ~isempty(varargin) && strcmpi(varargin{1}, "w_center")
    w_center = true;
else
    w_center = false;
end
    

% Utility function, same as wrapToPi
theta2pi = @(theta) theta - 2 * pi * floor((theta + pi) / (2 * pi));  % Wrap to [-pi, pi]

% Get motifs
motifs = cell(size(LC));
for i = 1:numel(LC)
    lc = LC{i};
    mtf = zeros(size(lc, 1), 4);
    [LCscore, LCPC, ~, LCmu] = myPCA(LC{i}, 2);  % (tPeriod, 2); (nParcels, 2); (1, nParcels)
    LCscore = squeeze(LCscore);  % Not necessary, just in case of any change in myPCA.
    LCmu = LCmu';  % (nParcels, 1)
    switch motif_type
        case "unscaled_centered_PC"
            mtf = [LCPC -LCPC];
        case "PC_max_projection"
            mtf = [LCPC .* max(LCscore), LCPC .* min(LCscore)] + LCmu;
        case "Point_max_projection"
            [~, max_idx] = max(LCscore);
            [~, min_idx] = min(LCscore);
            mtf = lc(:, [max_idx min_idx]);
        case "PC_max_cosine"
            theta = cart2pol(LCscore(:, 1), LCscore(:, 2));  % Angle when projecting to PC1-PC2 plane
            target_theta = [0, pi / 2, -pi, -pi / 2];   % Angle of the four extremes
            diff_theta = theta2pi(theta - target_theta);  % Angle to the targets (wrap to [-pi, pi])
            [~, min_idx] = min(abs(diff_theta));  % Find the one closest to the targets
            mtf(:, 1) = LCPC(:, 1) * LCscore(min_idx(1), 1) + LCmu;
            mtf(:, 2) = LCPC(:, 2) * LCscore(min_idx(2), 2) + LCmu;
            mtf(:, 3) = LCPC(:, 1) * LCscore(min_idx(3), 1) + LCmu;
            mtf(:, 4) = LCPC(:, 2) * LCscore(min_idx(4), 2) + LCmu;
        case "Point_max_cosine"
            theta = cart2pol(LCscore(:, 1), LCscore(:, 2));
            target_theta = [0, pi / 2, pi, -pi / 2];
            diff_theta = theta2pi(theta - target_theta);
            [~, min_idx] = min(abs(diff_theta));
            mtf = lc(:, min_idx);
        case "Slowest"
            v = vecnorm(lc(:, 2:end) - lc(:, 1:end - 1));  % (1, tPeriod - 1)
            [~, min_idx] = min(v);
            theta = cart2pol(LCscore(:, 1), LCscore(:, 2));  % Angle when projecting to PC1-PC2 plane
            target_theta = [0, pi / 2, -pi, -pi / 2] + theta(min_idx);   % Angle of the four extremes
            diff_theta = theta2pi(theta - target_theta);  % Angle to the targets (wrap to [-pi, pi])
            [~, min_idx] = min(abs(diff_theta));  % Find the one closest to the targets
            mtf = lc(:, min_idx);
    end
    if w_center
        motifs{i} = [mtf LCmu];
    else
        motifs{i} = mtf;
    end
end

end