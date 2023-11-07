function [] = step0()
    % 主函数
    [x, Fs] = audioread('.\audio_files\audiofile.wav'); % 读取音频文件
    
    % 计算STFT
    window = hamming(2048); % 定义窗口类型和长度
    nfft = 2048; % 定义FFT点数
    hop = length(window)/2; % 定义相邻窗的重叠样本数
    [S, f, t] = stft_analysis(x, window, nfft, hop, Fs); % STFT

    % 绘制时域信息
    figure;
    t = (0:length(x)-1)/Fs; % 创建时间向量
    plot(t, x);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Time Domain Signal');
    
    % 稍后绘制频域图
    
    % STFT 合成
    y = stft_synthesis(S, window, nfft, hop, Fs);

    % 绘制合成后的时域信号
    figure;
    t_y = (0:length(y)-1)/Fs; % 创建时间向量
    plot(t_y, y);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Synthesized Time Domain Signal');

    % 输出合成后的音频到文件
    outputFilename = '.\audio_files\synthesized_audio.wav';
    y=real(y);
    audiowrite(outputFilename, y, Fs); % 写入音频文件
end

function [S, f, t] = stft_analysis(x, window, nfft, hop, Fs)
    % 使用 MATLAB 内置的 stft 函数
    [S, f, t] = stft(x, Fs, 'Window', window, 'OverlapLength', hop, 'FFTLength', nfft);
end

function [x] = stft_synthesis(S, window, nfft, hop, Fs)
    % 使用 MATLAB 内置的 istft 函数
    x = istft(S, Fs, 'Window', window, 'OverlapLength', hop, 'FFTLength', nfft);
end

