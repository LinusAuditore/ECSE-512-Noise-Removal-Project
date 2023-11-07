function [] = step0()
    % 主函数
    [x, Fs] = audioread('.\audio_files\audiofile.wav'); % 读取音频文件
    
    % 计算STFT
    window = hamming(1024); % 定义窗口类型和长度
    nfft = 2048; % 定义FFT点数
    hop = 512; % 定义相邻窗的重叠样本数
    [S, f, t] = stft_analysis(x, window, nfft, hop, Fs); % STFT

    % 绘制时域信息
    figure;
    t = (0:length(x)-1)/Fs; % 创建时间向量
    plot(t, x);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Time Domain Signal');
    
    % 绘制STFT的大小（幅度谱）
    figure;
    S_mag = abs(S); % 计算幅度谱
    imagesc(t, f, 20*log10(S_mag)); % 绘制幅度谱的对数尺度
    axis xy; % 使y轴正向向下
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('Magnitude Spectrogram');
    
    % 绘制STFT的相位
    figure;
    S_phase = angle(S); % 计算相位谱
    imagesc(t, f, S_phase);
    axis xy;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('Phase Spectrogram');

    % 计算每个时间帧的幅度
    S_magnitude = abs(S);

    % 计算平均幅度（如果你想要显示所有时间帧的平均值）
    mean_magnitude = mean(S_magnitude, 2);

    % 绘制频率与幅度
    figure;
    plot(f, mean_magnitude);
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    title('Average Magnitude Spectrum');

    
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
    [S, f, t] = stft(x, Fs, 'Window', window, 'OverlapLength', length(window)-hop, 'FFTLength', nfft);
    % 转换 S 为单侧频谱
    S = S(1:nfft/2+1,:);
    % 调整频率向量以匹配单侧频谱
    f = f(1:nfft/2+1);
end


function [x] = stft_synthesis(S, window, nfft, hop, Fs)
    % 使用 MATLAB 内置的 istft 函数
    x = istft(S, Fs, 'Window', window, 'OverlapLength', length(window)-hop, 'FFTLength', nfft);
    full_length = length(x);  % 确定 x 的长度
    half_length = fix(full_length / 2);  % 计算一半长度的样本数，使用 fix 来获取整数值

    % 取前一半的时间的音频
    x = x(1:half_length);
end

