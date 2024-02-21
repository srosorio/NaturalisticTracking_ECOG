clear, clc

cd('F:\Matlab\IEEG\Scripts');

% add NSL toolbox to path
addpath('F:\Matlab\IEEG\Scripts\NSL');
addpath('F:\Matlab\fieldtrip-20210709');
ft_defaults;

% subject to use for dummy fieltrip structure
sub2plot = 'sub-02';    

% frequency band to analyze (SLB or HFB)
band2analyze = 'HFB'; 

% load audio file
[y,fs] = audioread('F:\Matlab\IEEG\Stimulus\audio.wav');
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
    
    % v corresponds to the signals in the 128 frequencies between 180 and 7142 Hz
    [music_v_all(:,:,idx),freqs]  = wav2aud2(musicsignal,[5 8 -2 0]);
    
    music_v  = music_v_all(:,:,idx); % detrend(music_v_all(:,:,idx));
    
    % get the envelope by averaging across frequencies
    envelope_music(idx,:)  = mean(music_v');
end

%%

% first, let's get the xcorr between the signal and copy of itself
one_env_music = envelope_music(1,:);
[r,lags]      = xcorr(one_env_music,one_env_music,400,'normalized');
figure(1),
plot(lags,r);

% now let's create the simulated signal
noise = wgn(1,length(one_env_music),1);
t     = linspace(0,30,length(one_env_music));     
   
if strcmpi(band2analyze,'SFB')
    sim_signal = one_env_music ./ max(abs(one_env_music)) + 1*noise;
else
%     sim_signal = one_env_music .* sin(85*2*pi*t) + noise;
    sim_signal = ( one_env_music ./ max(abs(one_env_music)) .* sin(85*2*pi*t) ) + 0.5*noise;
end

%% now load data from one subject and bring it down to epochs

addpath('F:\Matlab\fieldtrip-20210709');
ft_defaults;

% create a directory to save data
iEEG_dir = 'F:\Matlab\IEEG';
data_dir = [iEEG_dir,filesep,'Data'];

% get data from brainstorm database
cd(['F:\Matlab\brainstorm_db\iEEG\data\' sub2plot]);
FileTask = dir([sub2plot '*film*']);

% make sure to keep only good channels
load([cd,filesep,FileTask.name,filesep,'data_block001.mat'], 'ChannelFlag');
load([cd,filesep,FileTask.name,filesep,'channel.mat'], 'Channel');

% delete channels that won't be included from the ChannelFlag list
cd(['F:\Matlab\IEEG\' sub2plot]);

% load data
cfg                     = [];
cfg.dataset             = [sub2plot '.eeg'];
cfg.headerfile          = [sub2plot '.vmrk'];
cfg.trialfun            = 'ft_trialfun_general';
cfg.trialdef.eventtype  = 'Stimulus';
cfg.trialdef.eventvalue = 'speech';
cfg.trialdef.prestim    = 0;
cft.trialdef.poststim   = 30;
SpeechData              = ft_definetrial(cfg);

cfg.trialdef.eventtype  = 'Stimulus';
cfg.trialdef.eventvalue = 'music';
MusicData               = ft_definetrial(cfg);

%% define trials 
cfg = [];
cfg.dataset = [sub2plot '.eeg'];
cfg.trial           = 'all';
cfg.demean          = 'yes';
ThisData = ft_preprocessing(cfg);

% define trials according to music and speech segments. Create separate files
Speech_data = ft_redefinetrial(SpeechData,ThisData);
Music_data  = ft_redefinetrial(MusicData,ThisData);
%%
% resample data to match sampling rate in cochlear envelope
cfg = [];
cfg.resamplefs  = 200;
Speech_data     = ft_resampledata(cfg,Speech_data);
Music_data      = ft_resampledata(cfg,Music_data);

%% 
% create a dummy subject and replace actual data with the simulated signal.
% delete all other trials
for idx=1:length(Music_data.trial)
    if idx == 1
        Music_data.trial{idx} = repmat(sim_signal,length(Music_data.label),1);
        Music_data.time{idx}  = Music_data.time{idx}(1:length(sim_signal));
    else
        while length(Music_data.trial) > 1
            Music_data.trial(2) = [];
            Music_data.time(2)  = [];
        end
    end
end

%% 
% now apply ft filters to simulated datacfg = [];
cfg = [];
cfg.bpfilter        = 'yes';
cfg.bpfiltord       = 3;
if strcmpi(band2analyze,'SFB')
    cfg.bpfreq      = [1 8];
elseif strcmpi(band2analyze,'HFB')
    cfg.bpfreq      = [70 98];
    cfg.hilbert     = 'abs';
end

Music_data      = ft_preprocessing(cfg,Music_data);
%% 
% now get the cross-correlation between the processed simulated data and
% the envelope 

[r,lags]  = xcorr(zscore(one_env_music),zscore(Music_data.trial{1}(1,:)),400,'normalized');
figure(2)
plot(lags,r)