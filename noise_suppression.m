%% Main Function
function [] = noise_suppression()
    clc;
    clear;

    %% Obtain Clean Speech and Add Noise
    [input,Fs] = audioread('./audio_files/clean_speech.wav'); % get clean speech input and the sampling frequency Fs
    figure;
    plot(input);
    title('clean speech');
    % sound(input,Fs);
    % pause(10);
    
    % Add Noise
    snr = 10; %set Signal-Noise-Ratio to 10dB
    noised_signal = add_noise(input, "white noise", snr); % using white noise or pink noise
    figure;
    plot(noised_signal);
    title('clean speech+noise');
    % sound(noised_signal, Fs);
    % pause(10);
    
    % get the mean of the noise
    timeD_noise = noised_signal - input;
    timed_noise_mean = mean(timeD_noise);


    %% STFT
    window = 256;
    %effective analysis time interval is 38ms, window length is (time of speech)/38ms = 10s/38ms = 263
    noverlap = window/2;
    nfft = window;
    [stft_signal, f_signal, t_signal] = stft(noised_signal,Fs,'Window',hamming(window),'OverlapLength',noverlap,'FFTLength',nfft);
    [stft_noise, f_noise, t_noise] = stft(timeD_noise,Fs,'Window',hamming(window),'OverlapLength',noverlap,'FFTLength',nfft);


    %% Magnitude Averaging and Bias Removal
    magnitude_signal = abs(stft_signal); % get the magnitude of stft_signal
    
    % Magnitude Averaging
    magnitude_signal = magnitude_avg(magnitude_signal,3);
    stft_noise_mean = mean(abs(stft_noise)); % get [u1,u2,u3,...,uN];
    % Get H(e^jw) matrix
    Hejw = 1 - (stft_noise_mean ./ magnitude_signal); 
    
    % Half-Wave Rectification
    Hejw = 0.5 .* (Hejw + abs(Hejw));
    noise_removed_signal = Hejw .* stft_signal;
    
    % Reduce Noise Residual
    max_N_R = max(abs(stft_noise - stft_noise_mean .* (exp(1i * angle(stft_noise)))));
    noise_removed_signal = res_noise_reduction(noise_removed_signal, max_N_R);
    
    % Non-speech Activity Signal Attenuation
    avg_noise = stft_noise_mean .* exp(1i * angle(stft_signal));
    noise_removed_signal = non_speech_reduction(noise_removed_signal,avg_noise,stft_signal);
    
    
    %% ISTFT Function
    [output, Ts] = istft(noise_removed_signal,Fs,'Window',hamming(window),'OverlapLength',noverlap,'FFTLength',nfft);
    figure;
    output = real(output)./abs(real(output)) .* abs(output);
    % output = real(output);
    plot(output);
    title('Noise Removed Signal');
    sound(output,Fs);
    outputFilename = '.\audio_files\synthesized_audio.wav';
    audiowrite(outputFilename, output, Fs); % 写入音频文件
end

%% Add Noise
function [output] = add_noise(input,noise_type,snr)
    if noise_type == "white noise"
        output = awgn(input, snr,'measured');
    elseif noise_type == "pink noise"
        pn = pinknoise(size(input),'like', input);
        output = input + pn;
    else
        error('Unknown Noise Type');
    end

end


%% Magnitude Averaging
function [output] = magnitude_avg(signal_magnitude, avg_length)
    [~, length_column] = size(signal_magnitude);
    output = signal_magnitude;
    for column = 1:1:(length_column-avg_length+1)
        output(:,column) = mean(signal_magnitude(:,column:column+avg_length-1),2);
    end
end


%% Residual Noise Reduction
function [output] = res_noise_reduction(input, max_N_R)
    assert(size(input, 2) == length(max_N_R), 'max_N_R 大小不匹配');
    magnitude_input = abs(input);
    % output = magnitude_input;
    phase_input = angle(input);
    for column = 2:size(input, 2)-1
        max_N_R_column = max_N_R(1,column);
        for row = 1:size(input, 1)
            if magnitude_input(row,column) < max_N_R_column
                magnitude_input(row,column) = min(magnitude_input(row,column-1:column+1));
            end
        end
    end
    output = magnitude_input .* exp(1i * phase_input);
end


%% Non-Speech Signal Attenuation
function [output] = non_speech_reduction(input,avg_noise,original_signal)
    output = input;
    T = 20*log10((1/(2*pi))*sum(abs(input ./ avg_noise)));
    disp(T);
    for i=1:size(input,2)
        if T(1,i) < -12
            output(:,i) = 10^(-1.5) .* original_signal(:,i);
        end
    end
end