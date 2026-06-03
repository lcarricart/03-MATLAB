clc; clear;
mainTask1();

function mainTask1()

    omega = 2*pi*(0:0.01:100);   % rad/s

    % Compute FT for δ(t)
    ft_delta = computeFT_diracPulse(omega, 0);
    [mag_delta, phase_delta] = getMagPhaseSpectra_fromFTVec(ft_delta);

    % Compute FT for δ(t - 1 ms)
    ft_shift = computeFT_diracPulse(omega, 0.001);
    [mag_shift, phase_shift] = getMagPhaseSpectra_fromFTVec(ft_shift);

    % First figure: delta and shifted delta 
    figure;

    subplot(2,2,1);
    plot(omega, mag_delta);
    title('Magnitude spectrum δ(t)');

    subplot(2,2,2);
    plot(omega, phase_delta);
    title('Phase spectrum δ(t)');

    subplot(2,2,3);
    plot(omega, mag_shift);
    title('Magnitude spectrum δ(t - 1ms)');

    subplot(2,2,4);
    plot(omega, phase_shift);
    title('Phase spectrum δ(t - 1ms)');

    % Part (d): analytic vs linearity-based 

    % analytic spectra
    [mag_fx, phase_fx] = getMagPhaseSpectra_precalculated(omega);

    % linearity FT: δ(t) + δ(t - 1 ms)
    FT1 = computeFT_diracPulse(omega, 0);
    FT2 = computeFT_diracPulse(omega, 0.001);
    FT_combined = FT1 + FT2;

    [lin_mag, lin_phase] = getMagPhaseSpectra_fromFTVec(FT_combined);

    figure
    subplot(2,2,1);
    plot(omega, lin_mag);
    title('Combined Magnitude Spectrum');

    subplot(2,2,2);
    plot(omega, lin_phase);
    title('Combined Phase Spectrum');

    subplot(2,2,3);
    plot(omega, mag_fx);
    title('Precalculated Analytical Magnitude Spectrum');

    subplot(2,2,4);
    plot(omega, phase_fx);
    title('Precalculated Analytical Phase Spectrum');

    % ----- Part (e): New figure: f1(t) and f2(t) -----
    % f1(t) = δ(t) + δ(t − 1ms)
    F1_part1 = computeFT_diracPulse(omega, 0);
    F1_part2 = computeFT_diracPulse(omega, 0.001);
    F1_total = F1_part1 + F1_part2;
    [mag_f1, phase_f1] = getMagPhaseSpectra_fromFTVec(F1_total);

    % f2(t) = δ(t − 1ms) + δ(t − 2ms)
    F2_part1 = computeFT_diracPulse(omega, 0.001);
    F2_part2 = computeFT_diracPulse(omega, 0.002);
    F2_total = F2_part1 + F2_part2;
    [mag_f2, phase_f2] = getMagPhaseSpectra_fromFTVec(F2_total);

    figure
    subplot(2,2,1);
    plot(omega, mag_f1);
    title('Magnitude Spectrum of f1(t)');

    subplot(2,2,2);
    plot(omega, phase_f1);
    title('Phase Spectrum of f1(t)');

    subplot(2,2,3);
    plot(omega, mag_f2);
    title('Magnitude Spectrum of f2(t)');

    subplot(2,2,4);
    plot(omega, phase_f2);
    title('Phase Spectrum of f2(t)');

end

function ftVec = computeFT_diracPulse(omega, t0)
    ftVec = exp(-1j * omega * t0);
end

function [magSpec, phaseSpec] = getMagPhaseSpectra_fromFTVec(ftVec)
    magSpec = abs(ftVec);
    phaseSpec = atan2(imag(ftVec), real(ftVec));
end

function [magSpec, phaseSpec] = getMagPhaseSpectra_precalculated(omega)
    tau = 0.001;
    magSpec = 2 * cos(omega * tau / 2);
    phaseSpec = -omega * tau / 2;
end
