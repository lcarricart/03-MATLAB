function ss2_problem1_mine(n)
    % Variables definition
    T    = 0.0025;
    Ts   = T * n;
    fs   = 1 / Ts;

    % Time vector definition
    t           = 0:Ts:1;
    t_precise   = 0:0.00001:1; 

    % Signal definitions
    x           = 4*sin(2*pi*t) + cos(16*pi*t + pi/4);
    x_original  = 4*sin(2*pi*t_precise) + cos(16*pi*t_precise + pi/4);

    % Plotting (Remember that using plot() tries to unify the dots with straight lines. stem() is the right choice for digital signals)
    % Sampled Signal
    figure;
    subplot(2,1,1);
    hold on; grid on;
    stem(t, x, 'Color', 'r');
    stem(t_precise, x_original, ...
            'LineStyle','none', ...         % remove vertical lines
            'Marker','.', ...
            'MarkerFaceColor','none', ...   % no fill
            'Color','b');                   % red);       
    title("Sampled signal x[n] vs Original Signal x(t)");
    xlabel("n vs t"); ylabel("x[n] vs x(t)");


end