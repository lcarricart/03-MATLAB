function [L2D, s2D, tgrid] = path_length_2d_spline(T, p, grid_factor)
% PATH_LENGTH_2D_SPLINE  Line integral (arc length) of a smoothed 2D path.
% Uses cubic smoothing splines on x_meters and y_meters, then integrates
% speed = sqrt((dx/dt)^2 + (dy/dt)^2) over time/index.
%
% Inputs:
%   T           table with x_meters, y_meters and optionally a time column
%   p           smoothing parameter for csaps in [0,1]; p=1 interpolates
%               lower p applies more smoothing. Default 0.995
%   grid_factor evaluate derivatives on grid_factor*N samples. Default 10
%
% Outputs:
%   L2D   total 2D arc length [m]
%   s2D   cumulative arc length [m] defined on tgrid
%   tgrid evaluation grid corresponding to s2D
%
% Note: requires Curve Fitting Toolbox (csaps, fnder, fnval). If not
% available, falls back to polyline length.

    if nargin < 2 || isempty(p), p = 0.995; end
    if nargin < 3 || isempty(grid_factor), grid_factor = 10; end

    % Extract data
    x = T.x_meters(:);
    y = T.y_meters(:);

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
    good = isfinite(x) & isfinite(y) & isfinite(t);
    x = x(good); y = y(good); t = t(good);
    [t, iu] = unique(t, 'stable');
    x = x(iu); y = y(iu);

    if numel(t) < 2
        L2D = 0; s2D = 0; tgrid = t; return
    end

    % If csaps is missing, do polyline length and return
    if exist('csaps','file') ~= 2
        dxy = diff([x y], 1, 1);
        seg = hypot(dxy(:,1), dxy(:,2));
        s2D = [0; cumsum(seg)];
        L2D = s2D(end);
        tgrid = t;
        return
    end

    % fit smoothing splines and differentiate
    fx = csaps(t, x, p); dfx = fnder(fx);
    fy = csaps(t, y, p); dfy = fnder(fy);

    % integrate speed on a fine grid
    N = max(200, grid_factor * numel(t));
    tgrid = linspace(t(1), t(end), N);
    vx = fnval(dfx, tgrid);
    vy = fnval(dfy, tgrid);
    speed = hypot(vx, vy);

    s2D = cumtrapz(tgrid, speed);
    L2D = s2D(end);
end
