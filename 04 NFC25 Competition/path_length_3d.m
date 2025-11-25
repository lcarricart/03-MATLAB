function [L, s] = path_length_3d(T)
% PATH_LENGTH_3D  Arc length of a 3D polyline stored in table T
% Required columns: x_meters, y_meters, altitude_meters
%
% Outputs:
%   L  - total path length (meters)
%   s  - cumulative arc length at each point (meters), same length as rows in T (after cleaning)

    % Extract and shape
    x = T.x_meters(:);
    y = T.y_meters(:);
    z = T.altitude_meters(:);

    XYZ = [x y z];

    % Remove rows with NaNs or Infs
    good = all(isfinite(XYZ), 2);
    XYZ = XYZ(good, :);

    if size(XYZ,1) < 2
        L = 0;
        s = 0;
        return
    end

    % Remove consecutive duplicates to avoid zero-length segments
    keep = [true; any(diff(XYZ,1,1) ~= 0, 2)];
    XYZ = XYZ(keep, :);

    % Segment vectors and lengths
    dXYZ = diff(XYZ, 1, 1);
    seg = sqrt(sum(dXYZ.^2, 2));        % same as vecnorm(dXYZ,2,2)

    % Cumulative and total length
    s = [0; cumsum(seg)];
    L = s(end);
end
