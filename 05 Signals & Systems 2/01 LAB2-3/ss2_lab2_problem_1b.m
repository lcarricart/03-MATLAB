%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        %
%   SS2 Laboratory 2: Sampling Theorem - Problem 1       %
%                                                        %
%   Team Members: Georgii Molyboga  (2782258)            %
%                 Luciano Carricart (2782740)            %
%                 Mykyta  Kandyla   (2696614)            %
%                                                        %
%   Date: 07.05.2026                                     %
%                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a test signal of length 4
x = 0:1:4;

% Test MATLAB's fft
tic;
X_fft = fftshift(fft(x));
t_fft = toc;

% Test our myDFT
tic;
X_my = myDFT(x);
t_my = toc;

% Check maximum error to ensure correctness
err = max(abs(X_my - X_fft));

fprintf("my_DFT runtime: %.6f seconds\n", t_my);
fprintf("fft runtime:    %.6f seconds\n", t_fft);
fprintf("maximum error:  %.6e\n", err);

disp("X_my:");
disp(X_my);

disp("fftshift(fft(x)):");
disp(X_fft);