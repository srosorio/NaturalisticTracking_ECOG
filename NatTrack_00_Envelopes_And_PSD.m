clear, clc

cd('E:\Matlab\IEEG\Scripts');
% add NSL toolbox to path
addpath('E:\Matlab\IEEG\Scripts\NSL');

% create a directory to save data
iEEG_dir = 'E:\Matlab\IEEG';
data_dir = [iEEG_dir,filesep,'Data'];

% load audio file
[y,fs] = audioread('E:\Matlab\IEEG\Stimulus\audio.wav');
dt     = 1/fs;                   
tfull  = 0:dt:(length(y)*dt)-dt;  % time vector

% Cut full audio signal into speech and music segments
segstart  = 1;
seglength = 30;
totalsegments = floor((length(tfull)/48000) / 30);

% initiallize variables
speech_segments = cell(1,totalsegments);
music_segments  = cell(1,totalsegments);

for idx=1:totalsegments
    if mod(idx,2)==0
        % speech
        speech_segments{idx} = y(dsearchn(tfull',segstart):dsearchn(tfull',seglength*idx),:);
        segstart = segstart + seglength;
    elseif mod(idx,2)~=0
        % music 
        music_segments{idx}  = y(dsearchn(tfull',segstart):dsearchn(tfull',seglength*idx),:);
        segstart = segstart + seglength;
    end
end

% get rid of empty cells
music_segments(cellfun(@isempty,music_segments))   = [];
speech_segments(cellfun(@isempty,speech_segments)) = [];

% create a new time vector for segmented data
tseg = 0:dt:length(music_segments{1,1})*dt-dt; 

%% Get cochlear envelopes

for idx=1:length(music_segments)-1
    % resample signals
    musicsignal   = resample(music_segments{idx},tseg,16000);  %%%% resamples and then pads to 12 secs in total
    speechsignal  = resample(speech_segments{idx},tseg,16000);
    
    % v corresponds to the signals in the 128 frequencies between 180 and 7142 Hz
    [music_v,~]      = wav2aud2(musicsignal,[5 8 -2 0]);
    [speech_v,~]     = wav2aud2(speechsignal,[5 8 -2 0]);
    
    music_v  = detrend(music_v);
    speech_v = detrend(speech_v);
    
    % get the envelope by averaging across frequencies
    envelope_music(idx,:)  = mean(music_v');
    envelope_speech(idx,:) = mean(speech_v');
end

% save envelopes
save([data_dir,filesep,'envelopes_music.mat'],'envelope_music');
save([data_dir,filesep,'envelopes_speech.mat'],'envelope_speech');

%% Get PSDs

% low pass filter the data to look for slow frequency peaks
envelope_musicFILT  = lowpass(envelope_music',10,200,'Steepness',.95)';
envelope_speechFILT = lowpass(envelope_speech',10,200,'Steepness',.95)';

% get powers spectrum density for music
[psd_music, frex_music] = pwelch(envelope_musicFILT',[],[],[],250);
% normalize spectrum
psd_music = psd_music ./ sum(psd_music);

% plot
figure(1)
h1 = plot(frex_music,psd_music(:,1:6),'linew',1);
xlim([0 12]); ylim([0 .3]);
xlabel('Frequency (Hz)'); ylabel('Normalized power');
title('Music','FontWeight','normal');
set(gca,'FontSize',20);
box off, hold on

for idx=1:numel(h1)-1
    h1(idx).Color = [h1(idx).Color 0.3];
end

h2 = plot(frex_music,mean(psd_music(:,1:6),2),'linew',2); 
legend({'Segment 1', 'Segment 2', 'Segment 3', 'Segment 4', 'Segment 5', 'Segment 6','Average'},'FontSize',14,'Box','off');
% saveas(gca,'MusicStimuliPSD.tiff');

% now get the PSD for speech
[psd_speech, frex_speech] = pwelch(envelope_speechFILT',[],[],[],250);
% normalize spectrum
psd_speech = psd_speech ./ sum(psd_speech);

% plot
figure(2)
h1 = plot(frex_music,psd_speech,'linew',1);
xlim([0 12]); ylim([0 .3]);
xlabel('Frequency (Hz)'); ylabel('Normalized power');
title('Speech','FontWeight','normal');
set(gca,'FontSize',20);
box off, hold on

for idx=1:numel(h1)
    h1(idx).Color = [h1(idx).Color 0.3];
end

plot(frex_speech,mean(psd_speech,2),'linew',2,'Color', h2.Color);
legend({'Segment 1', 'Segment 2', 'Segment 3', 'Segment 4', 'Segment 5', 'Segment 6','Average'},'FontSize',14,'Box','off');
% saveas(gca,'SpeechStimuliPSD.tiff');

%% Plot raw acoustic signals
colors  = [0.7176 0.2745 1.0000; ...
           0.9412 0.5804 0.3373; ...
           0.4667 0.6745 0.1882];   
            
% plot music and speech signals

for idx=1:length(music_segments)-1
    figure(3)
    subplot(6,2,(idx*2)-1)
    ph = plot(tseg,music_segments{idx}(:,1),'color',colors(1,:));
    ylim([-1 1]); xlim([0 28.90]);
    box off
    if idx ~= 6
        set(gca,'XTick',[])
    end
    if idx == 1
        title('Signal','FontW','Normal')
    end    
    set(gca,'FontSize',20,'FontName','Arial')
end
xlabel('Time (s)')

for idx=1:length(speech_segments)
    figure(4)
    subplot(6,2,(idx*2)-1)
    ph = plot(tseg,speech_segments{idx}(:,1),'color',colors(2,:));
    ylim([-1 1]); xlim([0 28.90]);
    box off
    if idx ~= 6
        set(gca,'XTick',[])
    end
    if idx == 1
        title('Signal','FontW','Normal')
    end
    set(gca,'FontSize',20,'FontName','Arial')
end
xlabel('Time (s)')
%% 
% now plot the cochlear envelopes
for idx=1:length(music_segments)-1
    figure(3)
    subplot(6,2,idx*2)
    ph = plot(linspace(0,30,length(envelope_music)),envelope_music(idx,:),'color',colors(1,:));
    ylim([-1 1]); xlim([0 28.90]);
    box off
    if idx ~= 6
        set(gca,'XTick',[])
    end
    if idx == 1
        title('Envelope','FontW','Normal')
    end
    set(gca,'FontSize',20,'FontName','Arial')
end
xlabel('Time (s)')

% now plot the cochlear envelopes
for idx=1:length(speech_segments)
    figure(4)
    subplot(6,2,idx*2)
    ph = plot(linspace(0,30,length(envelope_speech)),envelope_speech(idx,:),'color',colors(2,:));
    ylim([-1 1]); xlim([0 28.90]);
    box off
    if idx ~= 6
        set(gca,'XTick',[])
    end
    if idx == 1
        title('Envelope','FontW','Normal')
    end
    set(gca,'FontSize',20,'FontName','Arial')
end
xlabel('Time (s)')

figure(3)
sh = suptitle('Music');
set(sh,'FontSize',24)
figure(4)
suptitle('Speech')
sh = suptitle('Speech');
set(sh,'FontSize',24)