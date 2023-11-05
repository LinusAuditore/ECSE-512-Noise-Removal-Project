function [] = step0()
    % 主函数
    [x, Fs] = audioread('.\audio_files\audiofile.wav'); % 读取音频文件

    % 绘制时域信息
    t = (0:length(x)-1)/Fs; % 创建时间向量
    plot(t, x);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Time Domain Signal');
    
    % 计算STFT
    [S, f, t] = stft_analysis(x, window, nfft, hop);
    
    % 绘制STFT的大小（幅度谱）
    S_mag = abs(S); % 计算幅度谱
    imagesc(t, f, 20*log10(S_mag)); % 绘制幅度谱的对数尺度
    axis xy; % 使y轴正向向下
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('Magnitude Spectrogram');
    colorbar; % 显示颜色条
    
    % 绘制STFT的相位
    S_phase = angle(S); % 计算相位谱
    imagesc(t, f, S_phase);
    axis xy;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('Phase Spectrogram');
    colorbar;
    
    % STFT 合成
    y = stft_synthesis(S, window, nfft, hop);

    % 绘制合成后的时域信号
    t_y = (0:length(y)-1)/Fs; % 创建时间向量
    plot(t_y, y);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Synthesized Time Domain Signal');

    % 输出合成后的音频到文件
    outputFilename = 'synthesized_audio.wav';
    audiowrite(outputFilename, y, Fs); % 写入音频文件
end

function [S, f, t] = stft_analysis(x, window, nfft, hop)
    % x: 输入信号
    % window: 窗函数
    % nfft: FFT的点数
    % hop: 相邻窗的重叠的样本数
    % S: STFT矩阵
    % f: 频率向量
    % t: 时间向量
    
    % 计算窗长度
    windowLen = length(window);
    
    % 确保信号是一个列向量
    x = x(:); 
    
    % 计算帧的数量
    numFrames = 1 + fix((length(x) - windowLen)/hop);
    
    % 初始化STFT矩阵
    S = complex(zeros(nfft, numFrames));
    
    % 对于每一帧进行FFT
    for k = 1:numFrames
        frame = x((k-1)*hop + (1:windowLen)).*window; % 应用窗函数
        S(:,k) = fft(frame, nfft); % 计算FFT
    end
    
    % 计算频率向量
    f = (0:nfft-1)*(Fs/nfft);
    
    % 计算时间向量
    t = (0:numFrames-1)*hop/Fs;
end

function [x] = stft_synthesis(S, window, nfft, hop)
    % S: STFT矩阵
    % window: 窗函数
    % nfft: FFT的点数
    % hop: 相邻窗的重叠的样本数
    % x: 合成后的信号
    
    % 计算窗长度
    windowLen = length(window);
    
    % 初始化合成信号
    xlen = windowLen + (size(S,2) - 1) * hop;
    x = zeros(xlen,1);
    
    % 初始化重叠-相加的窗口
    winAdd = zeros(xlen,1);
    
    % 对每一帧进行IFFT
    for k = 1:size(S,2)
        % 进行IFFT
        frame = ifft(S(:,k), nfft);
        
        % 取实部
        frame = real(frame(1:windowLen));
        
        % 重叠-相加
        x((k-1)*hop + (1:windowLen)) = x((k-1)*hop + (1:windowLen)) + frame.*window;
        winAdd((k-1)*hop + (1:windowLen)) = winAdd((k-1)*hop + (1:windowLen)) + window.^2;
    end
    
    % 处理重叠-相加的窗口，防止除以零
    winAdd(winAdd < 1e-6) = 1;
    
    % 正则化合成信号
    x = x ./ winAdd;
end
