% NaturalisticTracking_ECOG project
%
% This code estimates the crosscorrelation function between electrophysiological 
% signals and naturalistic acoustic (speech and music) signals. If specified, it
% also obtains a null distribution for the crosscorrelation function using
% a permutation procedure. 
%
% S.Osorio - 2023

clear, clc,

% create a directory to save data
iEEG_dir = 'E:\Matlab\IEEG';
data_dir = [iEEG_dir,filesep,'Data'];

% wheter to estimate observed xcorrelations (1 = yes, 0 = no)
observed_analysis    = 1;
% whether to estimate permutated xcorrelations (1 = yes, 0 = no)
permutation_analysis = 1;

band2analyze = 'HFB';  
sub2plot    = {'sub-02','sub-03','sub-05','sub-06','sub-10','sub-12','sub-16','sub-18', ...
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
load('E:\Matlab\IEEG\Scripts\envelopes_music.mat');
load('E:\Matlab\IEEG\Scripts\envelopes_speech.mat');

% whether to estimate crosscorrelation in small segments (1) or on the entire signal (0)
segestimation = 1;
fs            = 250;
n_trials      = length(AllDataStructuresFT{1,1}.trial);
n_conditions  = size(AllDataStructuresFT,2);                  
nperms        = factorial(n_trials);
maxlag        = 400;

%initialize cell arays for the data we need
[r_music,lag_music,r_speech,lag_speech] = deal(cell(1,length(sub2plot)));
[r_music_perm,lag_music_perm,r_speech_perm,lag_speech_perm] = deal(cell(1,length(sub2plot)));

% Observed croscorrelation
if observed_analysis == 1
    for sub_i=1:length(sub2plot)
        disp(['Estimating croscorrelation for subject ' num2str(sub_i)]);
        % number of electrodes (subject dependent)
        n_electrodes = size(AllDataStructuresFT{sub_i,1}.trial{1},1); % number of electrodes      
        % trim the EcoG data to match length of cochlear envelopes
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
                    % use sliding windows of 2 seconds
                    windowlength = 2 * fs;
                    kdx   = 1;
                    dtw_i = 1;
                    trialength = min([length(AllDataStructuresFT{sub_i,1}.trial{trial_i}(elec_i,:)) ...
                        length(envelope_speech(trial_i,:))]);  %length of current trial
                    while kdx < trialength
                        if trialength - kdx > windowlength
                            % speech
                            brain_signal     = AllDataStructuresFT{sub_i,1}.trial{trial_i}(elec_i,kdx:kdx+windowlength);
                            acoustic_signal  = envelope_speech(trial_i,kdx:kdx+windowlength);
                            % get all coefficients and lags (for z normalized signals)
                            [tempr,templags] = xcorr(zscore(brain_signal), ...
                                                     zscore(acoustic_signal),maxlag,'normalized');
                            % find max coefficient and its corresponding lag
                            r_speech{sub_i}(elec_i,trial_i,dtw_i)    = max(tempr);
                            lag_speech{sub_i}(elec_i,trial_i,dtw_i)  = templags(tempr == max(tempr));
                            % now music 
                            brain_signal     = AllDataStructuresFT{sub_i,2}.trial{trial_i}(elec_i,kdx:kdx+windowlength);
                            acoustic_signal  = envelope_music(trial_i,kdx:kdx+windowlength);
                            [tempr,templags] = xcorr(zscore(brain_signal), ...
                                                     zscore(acoustic_signal),maxlag,'normalized');
                            r_music{sub_i}(elec_i,trial_i,dtw_i)    = max(tempr);
                            lag_music{sub_i}(elec_i,trial_i,dtw_i)  = templags(tempr == max(tempr));
                            kdx = kdx+(windowlength/2);
                            dtw_i = dtw_i + 1;
                        else
                            % speech 
                            brain_signal     = AllDataStructuresFT{sub_i,1}.trial{trial_i}(elec_i,kdx:end);
                            acoustic_signal  = envelope_speech(trial_i,kdx:end);
                            [tempr,templags] = xcorr(zscore(brain_signal), ...
                                                     zscore(acoustic_signal),maxlag,'normalized');
                            r_speech{sub_i}(elec_i,trial_i,dtw_i)    = max(tempr);
                            lag_speech{sub_i}(elec_i,trial_i,dtw_i)  = templags(tempr == max(tempr));
                            % music 
                            brain_signal     = AllDataStructuresFT{sub_i,2}.trial{trial_i}(elec_i,kdx:end);
                            acoustic_signal  = envelope_music(trial_i,kdx:end);
                            [tempr,templags] = xcorr(zscore(brain_signal), ...
                                                     zscore(acoustic_signal),maxlag,'normalized');
                            r_music{sub_i}(elec_i,trial_i,dtw_i)    = max(tempr);
                            lag_music{sub_i}(elec_i,trial_i,dtw_i)  = templags(tempr == max(tempr));
                            kdx = trialength;
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
end

% Permuted Croscorrelations
if permutation_analysis == 1
    for sub_i=1:length(sub2plot)
        disp(['Estimating crosscorrelation for subject ' num2str(sub_i)]);
        n_electrodes = size(AllDataStructuresFT{sub_i,1}.trial{1},1); % number of electrodes
        % trim signal length if not already done
        if observed_analysis == 0
            for cond_i=1:n_conditions
                for trial_i=1:n_trials
                    AllDataStructuresFT{sub_i,cond_i}.trial{1,trial_i} = AllDataStructuresFT{sub_i,cond_i}.trial{1,trial_i}(:,1:length(envelope_speech));
                    AllDataStructuresFT{sub_i,cond_i}.time{trial_i}    = AllDataStructuresFT{sub_i,cond_i}.time{trial_i}(:,1:length(envelope_speech));
                end
            end
        end
        % permute n times
        for perm_i=1:nperms
            %disp(['Permutation ' num2str(perm_i)]);
            % shuffle trials for this permutation
            RandTrialOrder = randperm(numel(AllDataStructuresFT{sub_i,1}.trial));
            if segestimation == 1
                for trial_i=1:n_trials
                    for elec_i=1:n_electrodes
                        windowlength = 2 * fs; % use sliding windows of 2 seconds
                        kdx = 1; dtw_i = 1;
                        trialength = min([length(AllDataStructuresFT{sub_i,1}.trial{trial_i}(elec_i,:)) length(envelope_speech(trial_i,:))]);  %length of current trial               
                        while kdx < trialength
                            if trialength - kdx > windowlength
                                % speech
                                brain_signal     = AllDataStructuresFT{sub_i,1}.trial{trial_i}(elec_i,kdx:kdx+windowlength);
                                acoustic_signal  = envelope_speech(RandTrialOrder(trial_i),kdx:kdx+windowlength);
                                [tempr,templags] = xcorr(zscore(brain_signal), ... 
                                                         zscore(acoustic_signal),maxlag,'normalized');
                                r_speech_perm{sub_i}(elec_i,trial_i,dtw_i,perm_i) = max(tempr);
                                lag_speech_perm{sub_i}(elec_i,trial_i,dtw_i,perm_i)= templags(tempr == max(tempr));
                                % music
                                brain_signal     = AllDataStructuresFT{sub_i,2}.trial{trial_i}(elec_i,kdx:kdx+windowlength);
                                acoustic_signal  = envelope_music(RandTrialOrder(trial_i),kdx:kdx+windowlength);
                                [tempr,templags] = xcorr(zscore(brain_signal), ...
                                                         zscore(acoustic_signal),maxlag,'normalized');
                                r_music_perm{sub_i}(elec_i,trial_i,dtw_i,perm_i) = max(tempr);
                                lag_music_perm{sub_i}(elec_i,trial_i,dtw_i,perm_i)= templags(tempr == max(tempr));                         
                                kdx = kdx+(windowlength/2);
                                dtw_i = dtw_i + 1;
                            else
                                % speech
                                brain_signal     = AllDataStructuresFT{sub_i,1}.trial{trial_i}(elec_i,kdx:end);
                                acoustic_signal  = envelope_speech(RandTrialOrder(trial_i),kdx:end);
                                [tempr,templags] = xcorr(zscore(brain_signal), ...
                                                         zscore(acoustic_signal),maxlag,'normalized');
                                r_speech_perm{sub_i}(elec_i,trial_i,dtw_i,perm_i) = max(tempr);
                                lag_speech_perm{sub_i}(elec_i,trial_i,dtw_i,perm_i)= templags(tempr == max(tempr));
                                % and now music
                                brain_signal     = AllDataStructuresFT{sub_i,2}.trial{trial_i}(elec_i,kdx:end);
                                acoustic_signal  = envelope_music(RandTrialOrder(trial_i),kdx:end);
                                [tempr,templags] = xcorr(zscore(brain_signal), ...
                                                         zscore(acoustic_signal),maxlag,'normalized');
                                r_music_perm{sub_i}(elec_i,trial_i,dtw_i,perm_i) = max(tempr);
                                lag_music_perm{sub_i}(elec_i,trial_i,dtw_i,perm_i)= templags(tempr == max(tempr));
                                kdx = trialength;
                            end
                        end
                    end
                end
            else
                % if xcorr on entire segment
                for trial_i=1:n_trials
                    for elec_i=1:n_electrodes
                        % speech
                        brain_signal     = AllDataStructuresFT{sub_i,1}.trial{trial_i}(elec_i,:);
                        acoustic_signal  = envelope_speech(RandTrialOrder(trial_i),:);
                        [tempr,templags] = xcorr(zscore(brain_signal), ... 
                                                 zscore(acoustic_signal),maxlag,'normalized');
                        r_speech_perm{sub_i}(elec_i,trial_i,:,perm_i) = max(tempr);
                        lag_speech_perm{sub_i}(elec_i,trial_i,:,perm_i)= templags(tempr == max(tempr));
                        % music
                        brain_signal     = AllDataStructuresFT{sub_i,2}.trial{trial_i}(elec_i,:);
                        acoustic_signal  = envelope_music(RandTrialOrder(trial_i),:);
                        [tempr,templags] = xcorr(zscore(brain_signal), ...
                                                 zscore(acoustic_signal),maxlag,'normalized');
                        r_music_perm{sub_i}(elec_i,trial_i,:,perm_i) = max(tempr);
                        lag_music_perm{sub_i}(elec_i,trial_i,:,perm_i)= templags(tempr == max(tempr));
                    end
                end
            end
        end
    end
    
    if segestimation == 1
        disp(['saving xcorr_' band2analyze '_windowed_PERM.mat']);
        save([data_dir,filesep,'xcorr_' band2analyze '_windowed_PERM.mat'], ...
            'r_speech_perm','lag_speech_perm','r_music_perm', ...
            'lag_music_perm','band2analyze','sub2plot','AllChannelLabels','names4fields');
    else
        disp(['saving xcorr_' band2analyze '_PERM.mat']);
        save([data_dir,filesep,'xcorr_' band2analyze '_PERM.mat'], ...
            'r_speech_perm','lag_speech_perm','r_music_perm', ...
            'lag_music_perm','band2analyze','sub2plot','AllChannelLabels','names4fields');
    end
end
