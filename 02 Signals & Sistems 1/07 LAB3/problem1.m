%% Comments
% A Dirac pulse is not a filter signal, so my omega could be short-ranged (0 to 100 Hz)
% Y = fft(X) computes the discrete Fourier transform (DFT) of X using a fast Fourier transform (FFT) algorithm. Y is the same size as X.
% The vector omega is not for fft itself. It is the list of frequency points where you want to evaluate the (analytical) Fourier transform.

%% Start of the script
clc;
mainTask1();

%% Functions

% Exercise 1-B-i) Precompute the Fourier transform vector for a Dirac pulse ot time-shift t0, in the specific omega values specified in the vector
% Input: omega (a discrete vector of frequencies where we're interested in), t0 (time shift)
% Output: ftVec (Fourier Transform vector of a dirac with time shift t0)
function ftVec = computeFT_diracPulse(omega, t0) 
    % For a Dirac pulse, masking results in the following expression
    ftVec = exp(-1j.*omega .*t0);
end

% Exercise 1-B-ii) Decompose a FT into its magnitud and phase responses
% Input: ftVec (Fourier Transform of a function)
% Output: [magSpec, phaseSpec] (magnitude and phase response vectors)
function [magSpec, phaseSpec] = getMagPhaseSpectra_fromFTVec(ftVec)
    magSpec = abs(ftVec);
    
    % Applying SOHCAHTOA
    opposite = imag(ftVec);
    adjacent = real(ftVec);
    phaseSpec = atan2(opposite, adjacent);
end

% Exercise 1-B-iii) Compute the magnitude and phase of specific omega values, for the calculations done in exercise A)
% Input: omega (a discrete vector of frequencies where we're interested in)
% Output: [magSpec, phaseSpec] (magnitude and phase response vectors)
function [magSpec, phaseSpec] = getMagPhaseSpectra_precalculated(omega)
    t0 = 1e-3;
    FT = 1 + exp(-1j .* omega .* t0);
    magSpec = abs(FT);
    phaseSpec = angle(FT);
end

% Exercise 1-C) Plot the mag and phase responses of two Dirac pulses; a regular one, and a time-shifted one. Use subplot()
function mainTask1()
    %% Exercise 1-C)
    t0 = 1e-3;
    omega = 2*pi*(0:0.01:1000);
    fourierT1 = computeFT_diracPulse(omega, 0);
    fourierT2 = computeFT_diracPulse(omega, t0);
    [magSpec1, phaseSpec1] = getMagPhaseSpectra_fromFTVec(fourierT1);
    [magSpec2, phaseSpec2] = getMagPhaseSpectra_fromFTVec(fourierT2);
    
    figure(Name='Exercise 1) C)');
    subplot(2, 2, 1);
    stem(omega, magSpec1, 'LineStyle', 'none');
    title('Magnitude Response of the Regular Dirac Pulse');
    xlabel('Frequency (rad/s)');
    ylabel('Magnitude');
    xline(0);             % vertical line at ω = 0
    yline(0);             % horizontal line at 0 magnitude

    subplot(2, 2, 2);
    stem(omega, phaseSpec1, 'LineStyle', 'none');
    title('Phase Response of the Regular Dirac Pulse');
    xlabel('Angular Frequency (rad/s)');
    ylabel('Phase (radians)');
    xline(0);             
    yline(0);             

    subplot(2, 2, 3);
    stem(omega, magSpec2, 'LineStyle', 'none');
    title('Magnitude Response of the time-shifted Dirac Pulse at t0 = 1 ms');
    xlabel('Angular Frequency (rad/s)');
    ylabel('Magnitude');
    xline(0);             
    yline(0);

    subplot(2, 2, 4);
    stem(omega, phaseSpec2, 'LineStyle', 'none');
    title('Phase Response of the time-shifted Dirac Pulse at t0 = 1 ms');
    xlabel('Angular Frequency (rad/s)');
    ylabel('Phase (radians)');
    xline(0);             
    yline(0);             

    %% Exercise 1-D) I should see the cosine bending down.
    [ftPrecalc_magnitude, ftPrecalc_phase] = getMagPhaseSpectra_precalculated(omega);
    ftComputed1 = computeFT_diracPulse(omega, 0) + computeFT_diracPulse(omega , t0);
    [ftComputed_magnitude, ftComputed_phase] = getMagPhaseSpectra_fromFTVec(ftComputed1);

    figure(Name='Exercise 1) D)');
    subplot(2, 2, 1);
    stem(omega, ftPrecalc_magnitude, 'LineStyle', 'none');
    title('Magnitude Response of the Precalculated F(f(t))');
    xlabel('Angular Frequency (rad/s)');
    ylabel('Magnitude');
    xline(0);             
    yline(0);

    subplot(2, 2, 2);
    stem(omega, ftPrecalc_phase, 'LineStyle', 'none');
    title('Phase Response of the Precalculated F(f(t))');
    xlabel('Angular Frequency (rad/s)');
    ylabel('Phase (radians)');
    xline(0);             
    yline(0);

    subplot(2, 2, 3);
    stem(omega, ftComputed_magnitude, 'LineStyle', 'none');
    title('Magnitude Response of the Computed F(f(t))');
    xlabel('Angular Frequency (rad/s)');
    ylabel('Magnitude');
    xline(0);             
    yline(0);

    subplot(2, 2, 4);
    stem(omega, ftComputed_phase, 'LineStyle', 'none');
    title('Phase Response of the Computed F(f(t))');
    xlabel('Angular Frequency (rad/s)');
    ylabel('Phase (radians)');
    xline(0);             
    yline(0);

    fprintf("Precalculated phase and magnitue\n");
    fprintf('ftPrecalc_phase(98)      = %.6f rad\n', ftPrecalc_phase(98));
    fprintf('ftPrecalc_magnitude(340) = %.6f\n',     ftPrecalc_magnitude(340));
    
    fprintf("\nComputed phase and magnitue\n");
    fprintf('ftComputed_phase(98)      = %.6f rad\n', ftComputed_phase(98));
    fprintf('ftComputed_magnitude(340) = %.6f\n',     ftComputed_magnitude(340));

    %% Exercise 1-E)
    t1 = 2e-3;          % t0 was defined before
    ftComputed2 = computeFT_diracPulse(omega, t0) + computeFT_diracPulse(omega , t1);
    [ftComputed_magnitude2, ftComputed_phase2] = getMagPhaseSpectra_fromFTVec(ftComputed2);

    figure(Name='Exercise 1) E)');
    subplot(2, 2, 1);
    stem(omega, ftComputed_magnitude, 'LineStyle', 'none');
    title('Magnitude Response of f(t) = delta(t) + delta(t-1ms)');
    xlabel('Angular Frequency (rad/s)');
    ylabel('Magnitude');
    xline(0);             
    yline(0);

    subplot(2, 2, 2);
    stem(omega, ftComputed_phase, 'LineStyle', 'none');
    title('Phase Response of the f(t) = delta(t) + delta(t-1ms)');
    xlabel('Angular Frequency (rad/s)');
    ylabel('Phase (radians)');
    xline(0);             
    yline(0);

    subplot(2, 2, 3);
    stem(omega, ftComputed_magnitude2, 'LineStyle', 'none');
    title('Magnitude Response of f(t) = delta(t-1ms) + delta(t-2ms)');
    xlabel('Angular Frequency (rad/s)');
    ylabel('Magnitude');
    xline(0);             
    yline(0);

    subplot(2, 2, 4);
    stem(omega, unwrap(ftComputed_phase2), 'LineStyle', 'none'); % This unwrap() tells MATLAB that the angles are not jumping from quadrant to quadrant.
    title('Phase Response of f(t) = delta(t-1ms) + delta(t-2ms)');
    xlabel('Angular Frequency (rad/s)');
    ylabel('Phase (radians)');
    xline(0);             
    yline(0);
end