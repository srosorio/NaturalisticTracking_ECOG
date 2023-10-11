% NaturalisticTracking_ECOG project
%
% This code estimates the crosscorrelation function between electrophysiological 
% signals and naturalistic acoustic (speech and music) signals. 
%
% S.Osorio - 2023

clearvars -except time acoustic_signal brain_signal_HFB brain_signal_SFB
clc,

% create a directory to save data
iEEG_dir = 'F:\Matlab\IEEG';
cd(iEEG_dir)
data_dir = [iEEG_dir,filesep,'Data'];

band2analyze = 'SFB';  
sub2plot     = {'sub-02','sub-03','sub-05','sub-06','sub-10','sub-12','sub-16','sub-18', ...
                'sub-19','sub-20','sub-22','sub-24','sub-25','sub-26','sub-27','sub-34', ...
                'sub-36','sub-36HD','sub-39','sub-40','sub-45','sub-45HD','sub-46','sub-48', ...
                'sub-51','sub-54','sub-55','sub-58','sub-59','sub-60','sub-61','sub-63'};

% load neural data
if strcmpi(band2analyze,'SFB')
    load([data_dir,filesep,'fieldtrip_structures_SFB']);
elseif strcmpi(band2analyze,'HFB')
    load([data_dir,filesep,'fieldtrip_structures_HFB']);
end

% load acoustic envelopes
load('F:\Matlab\IEEG\Scripts\envelopes_music.mat');
load('F:\Matlab\IEEG\Scripts\envelopes_speech.mat');

segestimation = 0;      % whether to estimate xcorr in small overlapping segments (1) or on the entire signal (0)
fs            = 200;    % sample rate
n_trials      = length(AllDataStructuresFT{1,1}.trial);
n_conditions  = size(AllDataStructuresFT,2);
maxlag        = 400;    % maximum pos and neg lag for crosscorrelation

%initialize cell arays for the data we need
[r_music,lag_music,r_speech,lag_speech] = deal(cell(1,length(sub2plot)));

% Observed croscorrelation
for sub_i=1:length(sub2plot)
    disp(['Estimating croscorrelation for subject ' num2str(sub_i)]);
    n_electrodes = size(AllDataStructuresFT{sub_i,1}.trial{1},1);       % number of electrodes
    % trim end of EcoG signalk to match length of acoustic envelopes
    for cond_i=1:n_conditions
        for trial_i=1:n_trials
            AllDataStructuresFT{sub_i,cond_i}.trial{1,trial_i} = AllDataStructuresFT{sub_i,cond_i}.trial{1,trial_i}(:,1:length(envelope_speech));
            AllDataStructuresFT{sub_i,cond_i}.time{trial_i}    = AllDataStructuresFT{sub_i,cond_i}.time{trial_i}(:,1:length(envelope_speech));
        end
    end
    % croscorrelation if performed in sliding windows
    if segestimation == 1
        for trial_i=1:n_trials
            for elec_i=1:n_electrodes
                windowlength = 2 * fs; % 2-second window. Overlap is always windowlength / 2
                sample_i = 1;
                win_i    = 1;
                trialength = min([length(AllDataStructuresFT{sub_i,1}.trial{trial_i}(elec_i,:)) ...
                    length(envelope_speech(trial_i,:))]);  %length of current trial
                while sample_i < trialength
                    if trialength - sample_i > windowlength
                        % speech
                        brain_signal     = AllDataStructuresFT{sub_i,1}.trial{trial_i}(elec_i,sample_i:sample_i+windowlength);
                        acoustic_signal  = envelope_speech(trial_i,sample_i:sample_i+windowlength);
                        % get all coefficients and lags (for z normalized signals)
                        [tempr,templags] = xcorr(zscore(brain_signal), ...
                            zscore(acoustic_signal),maxlag,'normalized');
                        % find max coefficient and its corresponding lag
                        r_speech{sub_i}(elec_i,trial_i,win_i)    = max(tempr);
                        lag_speech{sub_i}(elec_i,trial_i,win_i)  = templags(tempr == max(tempr));
                        % now music
                        brain_signal     = AllDataStructuresFT{sub_i,2}.trial{trial_i}(elec_i,sample_i:sample_i+windowlength);
                        acoustic_signal  = envelope_music(trial_i,sample_i:sample_i+windowlength);
                        [tempr,templags] = xcorr(zscore(brain_signal), ...
                            zscore(acoustic_signal),maxlag,'normalized');
                        r_music{sub_i}(elec_i,trial_i,win_i)    = max(tempr);
                        lag_music{sub_i}(elec_i,trial_i,win_i)  = templags(tempr == max(tempr));
                        sample_i = sample_i+(windowlength/2);
                        win_i = win_i + 1;
                    else
                        % speech
                        brain_signal     = AllDataStructuresFT{sub_i,1}.trial{trial_i}(elec_i,sample_i:end);
                        acoustic_signal  = envelope_speech(trial_i,sample_i:end);
                        [tempr,templags] = xcorr(zscore(brain_signal), ...
                            zscore(acoustic_signal),maxlag,'normalized');
                        r_speech{sub_i}(elec_i,trial_i,win_i)    = max(tempr);
                        lag_speech{sub_i}(elec_i,trial_i,win_i)  = templags(tempr == max(tempr));
                        % music
                        brain_signal     = AllDataStructuresFT{sub_i,2}.trial{trial_i}(elec_i,sample_i:end);
                        acoustic_signal  = envelope_music(trial_i,sample_i:end);
                        [tempr,templags] = xcorr(zscore(brain_signal), ...
                            zscore(acoustic_signal),maxlag,'normalized');
                        r_music{sub_i}(elec_i,trial_i,win_i)    = max(tempr);
                        lag_music{sub_i}(elec_i,trial_i,win_i)  = templags(tempr == max(tempr));
                        sample_i = trialength;
                    end
                end
            end
        end
    % crosscorrelation if performed on the entire signal
    else
        for trial_i=1:n_trials
            for elec_i=1:n_electrodes
                % speech
                brain_signal     = AllDataStructuresFT{sub_i,1}.trial{trial_i}(elec_i,:);
                acoustic_signal  = envelope_speech(trial_i,:);
                [tempr,templags] = xcorr(zscore(brain_signal), ...
                    zscore(acoustic_signal),maxlag,'normalized');
                r_speech{sub_i}(elec_i,trial_i,:)    = max(tempr);
                lag_speech{sub_i}(elec_i,trial_i,:)  = templags(tempr == max(tempr));
                % music
                brain_signal     = AllDataStructuresFT{sub_i,2}.trial{trial_i}(elec_i,:);
                acoustic_signal  = envelope_music(trial_i,:);
                [tempr,templags] = xcorr(zscore(brain_signal), ...
                    zscore(acoustic_signal),maxlag,'normalized');
                r_music{sub_i}(elec_i,trial_i,:)    = max(tempr);
                lag_music{sub_i}(elec_i,trial_i,:)  = templags(tempr == max(tempr));
            end
        end
    end
end
if segestimation == 1
    save([data_dir,filesep,'xcorr_' band2analyze '_windowed.mat'], ...
        'r_speech','lag_speech','r_music','lag_music', ...
        'band2analyze','sub2plot','AllChannelLabels','names4fields');
else
    save([data_dir,filesep,'xcorr_' band2analyze '.mat'], ...
        'r_speech','lag_speech','r_music','lag_music', ...
        'band2analyze','sub2plot','AllChannelLabels','names4fields');
end