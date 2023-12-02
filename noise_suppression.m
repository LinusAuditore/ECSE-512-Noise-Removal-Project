clc;
clear;
%% Obtain Clean Speech and Add Noise
clean_speech = './clean speech/novel section.wav';
[input,Fs] = audioread(clean_speech); % get clean speech input and the sampling frequency Fs
% resampling to 16000Hz
[P,Q] = rat(16000/Fs);
abs(P/Q*Fs-16000);
Fs = 16000;
input = resample(input,P,Q);

N = length(input);
n = 0:N-1;
w = 2*n*pi/N;
% clean speech waveform plot
figure;
plot(input); 
xlabel('Samples','FontSize',15);
ylabel('Amplitude','FontSize',15);
title('Clean Speech Waveform Plot','FontSize',18);

% clean speech spectrogram
magnitude_clean = abs(fft(input));
len = round((length(magnitude_clean))/2);
figure;
plot(magnitude_clean(1:len));
xlabel('Frequency','FontSize',15);
ylabel('Magnitude','FontSize',15);
title('Clean Speech Spectrogram','FontSize',18);


sound(input,Fs);
pause(11);

% Add Noise
snr = 10; %set Signal-Noise-Ratio to 10dB
noised_signal = add_noise(input, "white noise", snr); % using white noise or pink noise
% noisy speech waveform plot
figure;
plot(noised_signal); % noisy speech waveform plot
xlabel('Samples','FontSize',15);
ylabel('Amplitude','FontSize',15);
title('Noisy Speech Waveform Plot','FontSize',18);

% noisy speech spectrogram
magnitude_noisy = abs(fft(noised_signal));
len = round((length(magnitude_noisy))/2);
figure;
plot(magnitude_noisy(1:len));
xlabel('Frequency','FontSize',15);
ylabel('Magnitude','FontSize',15);
title('Noisy Speech Spectrogram','FontSize',18);

sound(noised_signal, Fs);
pause(11);

% get the mean of the noise
timeD_noise = noised_signal - input;
timed_noise_mean = mean(timeD_noise);
%% STFT
window = 256;
%effective analysis time interval is 38ms, window length is (time of speech)/38ms = 10s/38ms = 263
noverlap = window/2;
nfft = window;
[stft_signal, f_signal, t_signal] = stft(noised_signal,Fs,'Window',hamming(window),'OverlapLength',noverlap,'FFTLength',nfft);
figure;
imagesc(t_signal, f_signal, 20*log10((abs(stft_signal(window/2:end,:)))));xlabel('Samples','FontSize',15); ylabel('Freqency','FontSize',15);
title('Spectrogram after STFT','FontSize',18);
colormap(gray);
colorbar;

% 3D waterfall spectrogram
% Short-time spectrum of noisy speech
figure; % figure 6
waterplot(stft_signal, f_signal, t_signal);

[stft_noise, f_noise, t_noise] = stft(timeD_noise,Fs,'Window',hamming(window),'OverlapLength',noverlap,'FFTLength',nfft);
%% just for comparison in report
% without magnitude averaging
magnitude_signal = abs(stft_signal); % get the magnitude of stft_signal
stft_noise_mean = mean(abs(stft_noise)); % get [u1,u2,u3,...,uN];

% Get H(e^jw) matrix
Hejw_without_magnitude_averaging = 1 - (stft_noise_mean ./ magnitude_signal);

% Half-Wave Rectification
Hejw_without_magnitude_averaging = 0.5 .* (Hejw_without_magnitude_averaging + abs(Hejw_without_magnitude_averaging));
noise_removed_signal_without_magnitude_averaging = Hejw_without_magnitude_averaging .* stft_signal;

% 3D waterfall spectrogram
% Short-time spectrum using bias removal and half-wave rectification
figure;% figure 7
waterplot(noise_removed_signal_without_magnitude_averaging, f_signal, t_signal);

%% Magnitude Averaging and Bias Removal
% Magnitude Averaging
mag_avg_num = 3;
magnitude_signal = magnitude_avg(magnitude_signal,mag_avg_num);

% 3D waterfall spectrogram
% Short-time spectrum of noisy-speech using three frame averaging
figure; % figure 8
waterplot(magnitude_signal, f_signal, t_signal);

% Get H(e^jw) matrix
Hejw = 1 - (stft_noise_mean ./ magnitude_signal);

% Half-Wave Rectification
Hejw = 0.5 .* (Hejw + abs(Hejw));
noise_removed_signal = Hejw .* stft_signal;

% 3D waterfall spectrogram
% Short-time spectrum using bias removal and half-wave rectification using three frame averaging
figure; % figure 9
waterplot(noise_removed_signal, f_signal, t_signal);

% Reduce Noise Residual
max_N_R = max(abs(stft_noise - stft_noise_mean .* (exp(1i * angle(stft_noise)))));
noise_removed_signal = res_noise_reduction(noise_removed_signal, max_N_R);

% Non-speech Activity Signal Attenuation
avg_noise = stft_noise_mean .* exp(1i * angle(stft_signal));
noise_removed_signal = non_speech_reduction(noise_removed_signal,avg_noise,stft_signal);

% spectrogram before ISTFT
figure;
imagesc(t_signal, f_signal, 20*log10((abs(noise_removed_signal(window/2:end,:)))));xlabel('Samples','FontSize',15); ylabel('Freqency','FontSize',15);
title('Spectrogram before ISTFT','FontSize',18);
colormap(gray);
colorbar;

% 3d waterfall spectrogram
% Short-time spectrum using bias removal, half-wave rectification, residual noise reduction, and nonspeech signal attenuation
figure; %figure 11
waterplot(noise_removed_signal, f_signal, t_signal);

%% ISTFT Function
[output, Ts] = istft(noise_removed_signal,Fs,'Window',hamming(window),'OverlapLength',noverlap,'FFTLength',nfft);

% noise removed speech waveform plot
figure;
% output = real(output)./abs(real(output)) .* abs(output);
output = real(output);
plot(output); % noise removed signal
xlabel('Samples','FontSize',15);
ylabel('Amplitude','FontSize',15);
title('Noise Removed Speech Waveform Plot','FontSize',18);

% noise removed speech spectrogram
magnitude_noise_removed = abs(fft(output));
len = round((length(magnitude_noise_removed))/2);
figure;
plot(magnitude_noise_removed(1:len));
xlabel('Frequency','FontSize',15);
ylabel('Magnitude','FontSize',15);
title('Noise Removed Speech Spectrogram','FontSize',18);
sound(output,Fs);
output_file_name = sprintf('./noise removed speech/output.wav', clean_speech, window, mag_avg_num);
audiowrite(output_file_name,output,Fs);
% FUNCTION
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
assert(size(input, 2) == length(max_N_R), 'max_N_R size not match');
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
for i=1:size(input,2)
    if T(1,i) < -12
        output(:,i) = 10^(-1.5) .* original_signal(:,i);
    end
end
end

function waterplot(s,f,t)
% Waterfall plot of spectrogram
    positive_freq_indices = f >= 0;
    f_positive = f(positive_freq_indices);
    s_positive = abs(s(positive_freq_indices,:)).^2;

    p = waterfall(f_positive, t, s_positive');
    % p = waterfall(f,t,abs(s)'.^2);
    set(gca,XDir="reverse",View=[-35 25]);
    xlabel("Frequency (Hz)", "FontSize",15);
    xlim([0, 4e3]);
    ylabel("Time (s)", "FontSize",15);
    ylim([0,10]); 
    zlabel("Amplitude", "FontSize",15);
    p.EdgeColor = [1,1,1];
    p.FaceColor = [0,0,0];
end