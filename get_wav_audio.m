%% Audio Recorder for subsequent processing
% Recording for 10 seconds
recObj = audiorecorder;
disp('Start Recording.');
recordblocking(recObj,10);
disp('Stop Recording');
% Play the recording
play(recObj);
% Get the recorded data
myrecording = getaudiodata(recObj);
% Plot the recorded data
plot(myrecording);
% Save the recorded data
filename = 'clean_speech.wav';
audiowrite(filename, myrecording, 8000); % sampling rate=8000

