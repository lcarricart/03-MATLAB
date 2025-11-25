%% Purpose of the program: Compute the Allan deviation across 10k seconds (2.7 hs) of data to understand what's the most effective window size that one should take for the averaging of a set. This project is directed towards improving the calibration of a sensor.

% % Load the data, keep only 10k values and the 6 relevant columns.
% T = readtable("Datalogger.xlsx", 'ReadVariableNames', true);
% T = T(1:10000, :);
% T = T(:, 1:7);
clc;

% Adjust the time array
time_ms = T.timestamp;
time_ms = time_ms - time_ms(1);                 % adjust the reference
time_s = time_ms / 1000;                        % time in seconds

%% Data Pre-Processing

% We need a uniform sampling for the Allan deviation
delta_t  = time_s(2) - time_s(1);        % seconds
t_start  = time_s(1);                    % first timestamp
t_end    = time_s(end);                  % last timestamp
t_uni    = t_start : delta_t : t_end;    % uniform grid in seconds

W = width(T);
T_uni = uniformData(t_uni, time_s, T, W);

%% Pre-Allan deviation

samplingPeriod = t_uni(2) - t_uni(1);           % The timestamps are separated approx every 2s
N              = numel(t_uni);                  % Number of samples
x_raw          = T_uni.gyroY;

% Choose m as powers of two up to the largest
maxPow = floor(log2(N/2));         
mList  = 2 .^ (0:maxPow);          % e.g. 1,2,4,...,4096
tau  = samplingPeriod * mList;                 % cluster times (s)

%% Perfom the Allan deviation

% Sweep over m values
for k = 1:numel(mList)
    m = mList(k);
    K = floor(N / m);              % number of complete clusters

    if K < 2                       
        sprintf('The analysed array of data does not allow Allan deviations')
        break
    end

    % Truncate so reshape works cleanly. (If N isnt an exact multiple of m, those trailing samples can't form a full cluster, so one simply drops them)
    x_trunc = x_raw(1:K*m);

    % Each column is one cluster of length m
    clusters = reshape(x_trunc, m, K);
    y        = mean(clusters, 1);           % cluster means y_k(τ)

    % Allan variance σ_A²(τ)  and deviation σ_A(τ)
    sigma2        = 0.5 * mean(diff(y).^2); % non-overlapping form
    adev(k)       = sqrt(sigma2);
end

%% Identify best averaging window
[adevMin, idxMin] = min(adev);      % deviation floor
tau0   = tau(idxMin);               % optimal τ (s)
m_opt  = mList(idxMin);             % optimal m  (samples)

%% ----------------- 5.  Plot & print results --------------------------
figure;
loglog(tau, adev, 'o-'); grid on, box on
xlabel('\tau  (s)'), ylabel('\sigma_A(\tau)  (sensor units)')
title('Allan Deviation of X-axis Accelerometer');
hold on
plot(tau0, adevMin, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
text(tau0*1.1, adevMin, sprintf('\\leftarrow  \\tau_0 = %.0f s', tau0));

fprintf('Optimal averaging window:\n');
fprintf('  τ₀  = %.0f s  (%.2f minutes)\n', tau0, tau0/60);
fprintf('  m₀  = %d samples\n', m_opt);
fprintf('  σ_A(τ₀) = %.4g (same units as input)\n', adevMin);