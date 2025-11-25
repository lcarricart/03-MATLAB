function [L2D, s2D] = path_length_2d(T)
% PATH_LENGTH_2D  Arc length of a 2D polyline using x_meters and y_meters
% Required columns in table T: x_meters, y_meters
%
% Outputs:
%   L2D - total 2D path length (meters)
%   s2D - cumulative arc length at each kept point (meters)

    % Extract and shape
    x = T.x_meters(:);
    y = T.y_meters(:);
    XY = [x y];

    % Remove rows with NaNs or Infs
    good = all(isfinite(XY), 2);
    XY = XY(good, :);

    % Must have at least two points
    if size(XY,1) < 2
        L2D = 0;
        s2D = 0;
        return
    end

    % Remove consecutive duplicates
    keep = [true; any(diff(XY,1,1) ~= 0, 2)];
    XY = XY(keep, :);

    % Segment vectors and lengths (use hypot for numeric stability)
    dXY = diff(XY, 1, 1);
    seg = hypot(dXY(:,1), dXY(:,2));

    % Cumulative and total length
    s2D = [0; cumsum(seg)];
    L2D = s2D(end);
end
