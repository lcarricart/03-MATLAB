function [L, s, tgrid] = path_length_3d_spline(T, p, grid_factor)
% PATH_LENGTH_3D_SPLINE  Line integral (arc length) along a smoothed 3D path.
% Inputs:
%   T           table with x_meters, y_meters, altitude_meters and optionally time
%   p           smoothing parameter for csaps in [0,1]; p=1 interpolates, lower smooths (default 0.995)
%   grid_factor evaluate derivatives on grid_factor*N points (default 10)
%
% Outputs:
%   L     total arc length [m]
%   s     cumulative arc length [m] on tgrid
%   tgrid evaluation grid for s

    if nargin < 2 || isempty(p), p = 0.995; end
    if nargin < 3 || isempty(grid_factor), grid_factor = 10; end

    x = T.x_meters(:); 
    y = T.y_meters(:); 
    z = T.altitude_meters(:);

    % choose parameter t: time if present, else sample index
    t = [];
    cand = {'timestamp','time','t','Time','Timestamp'};
    for k = 1:numel(cand)
        if ismember(cand{k}, T.Properties.VariableNames)
            t = T.(cand{k})(:);
            break
        end
    end
    if isempty(t), t = (0:numel(x)-1)'; end

    % remove bad rows and enforce strictly increasing t
    good = isfinite(x) & isfinite(y) & isfinite(z) & isfinite(t);
    x = x(good); y = y(good); z = z(good); t = t(good);
    [t, iu] = unique(t, 'stable');
    x = x(iu); y = y(iu); z = z(iu);
    if numel(t) < 2, L = 0; s = 0; tgrid = t; return; end

    % fit smoothing splines and differentiate
    fx = csaps(t, x, p); dfx = fnder(fx);
    fy = csaps(t, y, p); dfy = fnder(fy);
    fz = csaps(t, z, p); dfz = fnder(fz);

    % integrate speed = sqrt((dx/dt)^2+(dy/dt)^2+(dz/dt)^2)
    N = max(200, grid_factor * numel(t));
    tgrid = linspace(t(1), t(end), N);
    vx = fnval(dfx, tgrid);
    vy = fnval(dfy, tgrid);
    vz = fnval(dfz, tgrid);
    speed = sqrt(vx.^2 + vy.^2 + vz.^2);

    s = cumtrapz(tgrid, speed);
    L = s(end);
end
