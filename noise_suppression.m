clc;
clear;
%% Obtain Clean Speech and Add Noise
[input,Fs] = audioread('clean_speech.wav'); % get clean speech input and the sampling frequency Fs
% sound(input,Fs);
figure;
plot(input);
title('clean speech');
% Add Noise
snr = 10; %set Signal-Noise-Ratio to 10dB
noised_signal = add_noise(input, 'white noise', snr); % using white noise or pink noise
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
%% ISTFT Function
[output, Ts] = istft(noise_removed_signal,Fs,'Window',hamming(window),'OverlapLength',noverlap,'FFTLength',nfft);
figure;
plot(output);
title('Noise Removed Signal');
sound(output,Fs);
%% Add Noise
function [output]=add_noise(input,noise_type,snr)
if noise_type == 'white noise'
    output = awgn(input, snr,'measured');

elseif noise_type == 'pink noise'
    pn = pinknoise(size(input),'like',input);
    output = input + pn;
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

%% Bias Removal
function []=bias_removal(original_sig, noised_sig)

end
%% STFT with spectrogram function
% [s, f, t, p] = spectrogram(x, window, noverlap, nfft, Fs);
% figure;
% imagesc(t, f, 20*log10((abs(s))));xlabel('Samples'); ylabel('Freqency');
% title('使用spectrogram画出的短时傅里叶变换图形');
% colorbar;