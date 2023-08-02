% NaturalisticTracking_ECOG project
%
% This code estimates the crosscorrelation function between electrophysiological 
% signals and naturalistic acoustic (speech and music) signals. If specified, it
% also obtains a null distribution for the crosscorrelation function using
% a permutation procedure. 
%
% S.Osorio - 2023

clear, clc,
warning('off','all')
% initialize this session with a random seed, but keep track of it
rng('shuffle')
s = rng;

% create a directory to save data
iEEG_dir = '/autofs/cluster/transcend/sergio';
cd(iEEG_dir)
data_dir = [iEEG_dir,filesep,'Matlab_data'];

% wheter to estimate observed xcorrelations (1 = yes, 0 = no)
observed_analysis    = 0;
% whether to estimate permutated xcorrelations (1 = yes, 0 = no)
permutation_analysis = 1;
% frequency band to analyze
band2analyze = 'SFB';  
subjects     = {'sub-02','sub-03','sub-05','sub-06','sub-10','sub-12','sub-16','sub-18', ...
                'sub-19','sub-20','sub-22','sub-24','sub-25','sub-26','sub-27','sub-34', ...
                'sub-36','sub-36HD','sub-39','sub-40','sub-45','sub-45HD','sub-46','sub-48', ...
                'sub-51','sub-54','sub-55','sub-58','sub-59','sub-60','sub-61','sub-63'};
% run permutations only for this subjct
sub2permute  = 'sub-02';
sub2perm_idx = strcmpi(sub2permute,subjects);

% load neural data
if strcmpi(band2analyze,'SFB')
    load([data_dir,filesep,'fieldtrip_structures_SFB']);
elseif strcmpi(band2analyze,'HFB')
    load([data_dir,filesep,'fieldtrip_structures_HFB']);
end

% load acoustic envelopes
load([data_dir,filesep,'envelopes_music.mat']);
load([data_dir,filesep,'envelopes_speech.mat']);

% whether to estimate crosscorrelation in small segments (1) or on the entire signal (0)
segestimation = 0;
perm_type     = 'wn'; % 'wn' white noise, 'ts' trial shuffling
fs            = 250;
n_trials      = length(AllDataStructuresFT{1,1}.trial);
n_conditions  = size(AllDataStructuresFT,2);
maxlag        = 400;

% number of permutations depends on type of permutation
if strcmpi(perm_type,'ts')
    nperms = factorial(n_trials);
elseif strcmpi(perm_type,'wn')
    nperms = 1000;
end

% set filters according to the frequency band of interest
if strcmpi(band2analyze,'SFB')
    bands4filt = [1,8];
elseif strcmpi(band2analyze,'HFB')
    bands4filt = [70,120];
end 

%initialize cell arays for the data we need
[r_music_perm,lag_music_perm,r_speech_perm,lag_speech_perm] = deal(cell(1,1));
% get only th esubject of interest from the data structure
this_subject = {AllDataStructuresFT{sub2perm_idx,1},AllDataStructuresFT{sub2perm_idx,2}};
n_electrodes = size(this_subject{1,1}.trial{1},1); % number of electrodes

tic()
% All set, let's get our dummy croscorrelation values and lags
if permutation_analysis == 1
        disp(['NatTrack >>> Estimating null distribution for crosscorrelation coefficients for subject ' sub2permute]); 
        % trim signal length if not already done
        if observed_analysis == 0
            for cond_i=1:n_conditions
                for trial_i=1:n_trials
                    this_subject{1,cond_i}.trial{1,trial_i} = this_subject{1,cond_i}.trial{1,trial_i}(:,1:length(envelope_speech));
                    this_subject{1,cond_i}.time{trial_i}    = this_subject{1,cond_i}.time{trial_i}(:,1:length(envelope_speech));
                end
            end
        end
        if strcmpi(perm_type,'ts')
            for perm_i=1:nperms
                % shuffle trials for this permutation
                RandTrialOrder = randperm(numel(this_subject{1,1}.trial));
                % estimates permuations by randomly shuffling trial labels and using overlaping sliding windows
                if segestimation == 1
                    for trial_i=1:n_trials
                        for elec_i=1:n_electrodes
                            windowlength = 2 * fs; % use sliding windows of 2 seconds
                            kdx = 1; win_i = 1;
                            trialength = min([length(this_subject{1,1}.trial{trial_i}(elec_i,:)) length(envelope_speech(trial_i,:))]);  %length of current trial
                            while kdx < trialength
                                if trialength - kdx > windowlength
                                    % speech
                                    brain_signal     = this_subject{1,1}.trial{trial_i}(elec_i,kdx:kdx+windowlength);
                                    acoustic_signal  = envelope_speech(RandTrialOrder(trial_i),kdx:kdx+windowlength);
                                    [tempr,templags] = xcorr(zscore(brain_signal), ...
                                        zscore(acoustic_signal),maxlag,'normalized');
                                    r_speech_perm{1}(elec_i,trial_i,win_i,perm_i) = max(tempr);
                                    lag_speech_perm{1}(elec_i,trial_i,win_i,perm_i)= templags(tempr == max(tempr));
                                    % music
                                    brain_signal     = this_subject{1,2}.trial{trial_i}(elec_i,kdx:kdx+windowlength);
                                    acoustic_signal  = envelope_music(RandTrialOrder(trial_i),kdx:kdx+windowlength);
                                    [tempr,templags] = xcorr(zscore(brain_signal), ...
                                        zscore(acoustic_signal),maxlag,'normalized');
                                    r_music_perm{1}(elec_i,trial_i,win_i,perm_i) = max(tempr);
                                    lag_music_perm{1}(elec_i,trial_i,win_i,perm_i)= templags(tempr == max(tempr));
                                    kdx = kdx+(windowlength/2);
                                    win_i = win_i + 1;
                                else
                                    % speech
                                    brain_signal     = this_subject{1,1}.trial{trial_i}(elec_i,kdx:end);
                                    acoustic_signal  = envelope_speech(RandTrialOrder(trial_i),kdx:end);
                                    [tempr,templags] = xcorr(zscore(brain_signal), ...
                                        zscore(acoustic_signal),maxlag,'normalized');
                                    r_speech_perm{1}(elec_i,trial_i,win_i,perm_i) = max(tempr);
                                    lag_speech_perm{1}(elec_i,trial_i,win_i,perm_i)= templags(tempr == max(tempr));
                                    % and now music
                                    brain_signal     = this_subject{1,2}.trial{trial_i}(elec_i,kdx:end);
                                    acoustic_signal  = envelope_music(RandTrialOrder(trial_i),kdx:end);
                                    [tempr,templags] = xcorr(zscore(brain_signal), ...
                                        zscore(acoustic_signal),maxlag,'normalized');
                                    r_music_perm{1}(elec_i,trial_i,win_i,perm_i) = max(tempr);
                                    lag_music_perm{1}(elec_i,trial_i,win_i,perm_i)= templags(tempr == max(tempr));
                                    kdx = trialength;
                                end
                            end
                        end
                    end
                % otherwise we estimate permutations by label shuffling but using the full 30 second segment
                else
                    for trial_i=1:n_trials
                        for elec_i=1:n_electrodes
                            if strcmpi(perm_type,'ts')
                                % speech
                                brain_signal     = this_subject{1,1}.trial{trial_i}(elec_i,:);
                                acoustic_signal  = envelope_speech(RandTrialOrder(trial_i),:);
                                [tempr,templags] = xcorr(zscore(brain_signal), ...
                                    zscore(acoustic_signal),maxlag,'normalized');
                                r_speech_perm{1}(elec_i,trial_i,:,perm_i) = max(tempr);
                                lag_speech_perm{1}(elec_i,trial_i,:,perm_i)= templags(tempr == max(tempr));
                                % music
                                brain_signal     = this_subject{1,2}.trial{trial_i}(elec_i,:);
                                acoustic_signal  = envelope_music(RandTrialOrder(trial_i),:);
                                [tempr,templags] = xcorr(zscore(brain_signal), ...
                                    zscore(acoustic_signal),maxlag,'normalized');
                                r_music_perm{1}(elec_i,trial_i,:,perm_i) = max(tempr);
                                lag_music_perm{1}(elec_i,trial_i,:,perm_i)= templags(tempr == max(tempr));
                            end
                        end
                    end
                end
            end
        elseif strcmpi(perm_type,'wn')
            tic()
            % do permutations
            disp('NatTrack >>> Initializing permutations using white noise')
            for perm_i=1:nperms
                % get the time vector
                time_vector  = this_subject{1,cond_i}.time{trial_i};
                % create 6 white noise segments of 30 seconds per condition
                wn_musicsignal  = resample(wgn(6,length(envelope_speech),1)',time_vector,16000)';
                wn_speechsignal = resample(wgn(6,length(envelope_speech),1)',time_vector,16000)';
                % now we apply the cochlear filter to the wn signals as we did the real signals
                for sig_i=1:size(wn_speechsignal,1)
                    % v corresponds to the signals in the 128 frequencies between 180 and 7142 Hz
                    [music_v,~]  = wav2aud2(wn_musicsignal(sig_i,:),[5 8 -2 0]);
                    [speech_v,~] = wav2aud2(wn_speechsignal(sig_i,:),[5 8 -2 0]);
                    % detrend the signals
                    music_v  = detrend(music_v);
                    speech_v = detrend(speech_v);
                    % get the envelope by averaging across frequencies
                    envelope_music(sig_i,:)  = mean(music_v');
                    envelope_speech(sig_i,:) = mean(speech_v');
                end
                % now apply preprocessing filters to the fake acoustic envelopes
                if strcmpi(band2analyze,'SFB')
                    envelope_music  = bandpass(envelope_music',bands4filt,200)';
                    envelope_speech = bandpass(envelope_speech',bands4filt,200)';
                elseif strcmpi(band2analyze,'HFB')
                    envelope_music  = abs(hilbert(bandpass(envelope_music',bands4filt,200)))';
                    envelope_speech = abs(hilbert(bandpass(envelope_speech',bands4filt,200)))';
                end
                for trial_i=1:n_trials
                    for elec_i=1:n_electrodes
                        % speech
                        acoustic_signal  = envelope_speech(trial_i,:);
                        brain_signal     = this_subject{1,1}.trial{trial_i}(elec_i,:);
                        [tempr,templags] = xcorr(zscore(brain_signal), ...
                            zscore(acoustic_signal),maxlag,'normalized');
                        r_speech_perm{1}(elec_i,trial_i,perm_i) = max(tempr);
                        lag_speech_perm{1}(elec_i,trial_i,perm_i)= templags(tempr == max(tempr));
                        % music
                        acoustic_signal  = envelope_music(trial_i,:);
                        brain_signal     = this_subject{1,2}.trial{trial_i}(elec_i,:);
                        [tempr,templags] = xcorr(zscore(brain_signal), ...
                            zscore(acoustic_signal),maxlag,'normalized');
                        r_music_perm{1}(elec_i,trial_i,perm_i) = max(tempr);
                        lag_music_perm{1}(elec_i,trial_i,perm_i)= templags(tempr == max(tempr));
                    end
                end
            end
        end

   toc()
   % save data
   disp(['NatTrack >>> Done! Saving files for subject ' sub2permute])
   if segestimation == 1
        disp(['saving xcorr_' band2analyze '_windowed_PERM.mat']);
        save([data_dir,filesep,subjects{sub2perm_idx} '_xcorr_' band2analyze '_windowed_PERM.mat'], ...
            'r_speech_perm','lag_speech_perm','r_music_perm', ...
            'lag_music_perm','band2analyze','sub2plot','AllChannelLabels','names4fields','s');
    else
        if strcmpi(perm_type,'ts')
            disp(['saving xcorr_' band2analyze '_trialshuffle_PERM.mat']);
            save([data_dir,filesep,subjects{sub2perm_idx} '_xcorr_' band2analyze '_trialshuffle_PERM.mat'], ...
                'r_speech_perm','lag_speech_perm','r_music_perm', ...
                'lag_music_perm','band2analyze','sub2plot','AllChannelLabels','names4fields','s');
        else
            disp(['saving xcorr_' band2analyze '_whitenoise_PERM.mat']);
            save([data_dir,filesep,subjects{sub2perm_idx} '_xcorr_' band2analyze '_whitenoise_PERM.mat'], ...
                'r_speech_perm','lag_speech_perm','r_music_perm', ...
                'lag_music_perm','band2analyze','subjects','AllChannelLabels','names4fields','s');     
        end
    end
end