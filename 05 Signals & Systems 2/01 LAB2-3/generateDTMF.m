function y = generateDTMF(digits)
    fs = 8000;              % 8 kHz sample rate
    tone_duration = 0.075;  % 75 ms
    break_duration = 0.030; % 30 ms
    
    % Create time vectors
    t_tone = 0:(1/fs):(tone_duration - 1/fs);
    silence = zeros(1, round(break_duration * fs));
    
    y = []; % Initialize empty output array
    
    % Loop through each character in the input string
    for i = 1:length(digits)
        digit = digits(i);
        
        % Determine frequencies based on the DTMF table
        switch digit
            case '1', f1 = 697; f2 = 1209;
            case '2', f1 = 697; f2 = 1336;
            case '3', f1 = 697; f2 = 1477;
            case '4', f1 = 770; f2 = 1209;
            case '5', f1 = 770; f2 = 1336;
            case '6', f1 = 770; f2 = 1477;
            case '7', f1 = 852; f2 = 1209;
            case '8', f1 = 852; f2 = 1336;
            case '9', f1 = 852; f2 = 1477;
            case '*', f1 = 941; f2 = 1209;
            case '0', f1 = 941; f2 = 1336;
            case '#', f1 = 941; f2 = 1477;
            otherwise
                error('Invalid DTMF digit. Use 0-9, *, or #.');
        end
        
        % Generate the dual-tone signal for the current digit
        tone = sin(2*pi*f1*t_tone) + sin(2*pi*f2*t_tone);
        
        % Append the tone and the subsequent silence break to the total signal
        y = [y, tone, silence]; 
    end
end