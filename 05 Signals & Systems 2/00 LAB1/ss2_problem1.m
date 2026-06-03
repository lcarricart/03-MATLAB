%% SS2 Laboratory 1: Sampling Theorem

function ss2_problem1(n)
    % Variables definition
    T = 0.0025;
    t = 0:T:1;
    Ts = T * n;

    % Target function definition
    x = 4*sin(2*pi*t) + cos(pi/4 + 16*pi*t); 

    % Impose identical scaling
    x_sampled = downsample(x, n);
    t_sampled = 0:Ts:(length(x_sampled)-1)*Ts;

    % Mathematical signal reconstruction
    x_reconstructed = zeros(size(t));

    % sinc Interpolation
    for k = 0:length(x_sampled)-1
        x_current = x_sampled(k+1) * sinc((t - k*Ts) / Ts);
        x_reconstructed = x_reconstructed + x_current;
    end
    
    % Plots
    figure;
    subplot(2,1,1);
    plot(t, x);
    hold on;
    title(['Original Signal + Samples (n = ', num2str(n), ')']);
    ylabel("Amplitude"); xlabel("Time (s)");

    stem(t_sampled, x_sampled)

    subplot(2,1,2);
    plot(t, x_reconstructed);
    title("Reconstructed Signal");
    ylabel("Amplitude"); xlabel("Time (s)");
end