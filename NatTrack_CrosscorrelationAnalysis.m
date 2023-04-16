clear, clc

% create a directory to save data
iEEG_dir = 'E:\Matlab\IEEG';
data_dir = [iEEG_dir,filesep,'Data'];

band2analyze = 'HFB';  
sub2plot    = {'sub-02','sub-03','sub-05','sub-06','sub-10','sub-12','sub-16','sub-18', ...
                'sub-19','sub-20','sub-22','sub-24','sub-25','sub-26','sub-27','sub-34', ...
                'sub-36','sub-36HD','sub-39','sub-40','sub-45','sub-45HD','sub-46','sub-48', ...
                'sub-51','sub-54','sub-55','sub-58','sub-59','sub-60','sub-61','sub-63'};

% whether to estimate crosscorrelation in small segments (1) or on the entire signal (0)
segestimation = 1;
fs            = 250;
nperms        = 100;
maxlag        = 400;

% load neural data
if strcmpi(band2analyze,'SFB')
    load([data_dir,filesep,'fieldtrip_structures_SFB']);
elseif strcmpi(band2analyze,'HFB')
    load([data_dir,filesep,'fieldtrip_structures_HFB']);
end

% load acoustic envelopes
load('E:\Matlab\IEEG\Scripts\envelopes_music.mat');
load('E:\Matlab\IEEG\Scripts\envelopes_speech.mat');

% Observed croscorrelation
for sub_i=1:length(sub2plot)
    disp(['Estimating croscorrelation for subject ' num2str(sub_i)]);
     
    n_trials     = length(AllDataStructuresFT{sub_i,1}.trial);    % number of trials
    n_electrodes = size(AllDataStructuresFT{sub_i,1}.trial{1},1); % number of electrodes
    
    % SPEECH (first column of AllDatastructuresFT)
    % trim the EcoG data to match length of cochlear envelopes
    for idx=1:n_trials
        AllDataStructuresFT{sub_i,1}.trial{1,idx} = AllDataStructuresFT{sub_i,1}.trial{1,idx}(:,1:length(envelope_speech));
        AllDataStructuresFT{sub_i,1}.time{idx}    = AllDataStructuresFT{sub_i,1}.time{idx}(:,1:length(envelope_speech));
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
                        % z normalize the data 
                        brain_signal     = AllDataStructuresFT{sub_i,1}.trial{trial_i}(elec_i,kdx:kdx+windowlength);
                        acoustic_signal  = envelope_speech(trial_i,kdx:kdx+windowlength);
                        % get all coefficients and lags (for z normalized signals)
                        [tempr,templags] = xcorr(zscore(brain_signal), ...
                                                 zscore(acoustic_signal),maxlag,'normalized');
                        % find max coefficient and its corresponding lag
                        rdata_speech{sub_i}(elec_i,trial_i,dtw_i)    = max(tempr);
                        lagdata_speech{sub_i}(elec_i,trial_i,dtw_i)  = templags(find(tempr == max(tempr)));                                        
                        kdx = kdx+(windowlength/2);
                        dtw_i = dtw_i + 1;
                    else       
                        brain_signal     = AllDataStructuresFT{sub_i,1}.trial{trial_i}(elec_i,kdx:end);
                        acoustic_signal  = envelope_speech(trial_i,kdx:end);
                        [tempr,templags] = xcorr(zscore(brain_signal), ...
                                                 zscore(acoustic_signal),maxlag,'normalized');
                        rdata_speech{sub_i}(elec_i,trial_i,dtw_i)    = max(tempr);
                        lagdata_speech{sub_i}(elec_i,trial_i,dtw_i)  = templags(find(tempr == max(tempr)));
                        kdx = trialength;
                    end
                end
            end
        end  
    % crosscorrelation if performed on the entire signal
    else    
        for trial_i=1:n_trials
            for elec_i=1:n_electrodes
                brain_signal     = AllDataStructuresFT{sub_i,1}.trial{trial_i}(elec_i,:);
                acoustic_signal  = envelope_speech(trial_i,:);
                [tempr,templags] = xcorr(zscore(brain_signal), ...
                                         zscore(acoustic_signal),maxlag,'normalized');
                rdata_speech{sub_i}(elec_i,trial_i,:)    = max(tempr);
                lagdata_speech{sub_i}(elec_i,trial_i,:)  = templags(find(tempr == max(tempr)));
            end
        end
    end
    
    % MUSIC ( all the same, but for the second column of AllDatastructuresFT)
    for idx=1:length(AllDataStructuresFT{sub_i,2}.trial)
        AllDataStructuresFT{sub_i,2}.trial{1,idx} = AllDataStructuresFT{sub_i,2}.trial{1,idx}(:,1:length(envelope_music));
        AllDataStructuresFT{sub_i,2}.time{idx}    = AllDataStructuresFT{sub_i,2}.time{idx}(:,1:length(envelope_speech));
    end 
    
    if segestimation == 1
        for trial_i=1:n_trials
            for elec_i=1:n_electrodes
                windowlength = 2 * fs; % use sliding windows of 2 seconds
                kdx   = 1; 
                dtw_i = 1;  
                trialength = min([length(AllDataStructuresFT{sub_i,1}.trial{trial_i}(elec_i,:)) ...
                                  length(envelope_speech(trial_i,:))]);  %length of current trial
                while kdx < trialength
                    if trialength - kdx > windowlength
                        brain_signal     = AllDataStructuresFT{sub_i,2}.trial{trial_i}(elec_i,kdx:kdx+windowlength);
                        acoustic_signal  = envelope_music(trial_i,kdx:kdx+windowlength);
                        [tempr,templags] = xcorr(zscore(brain_signal), ... 
                                                 zscore(acoustic_signal),maxlag,'normalized');
                        rdata_music{sub_i}(elec_i,trial_i,dtw_i)    = max(tempr);
                        lagdata_music{sub_i}(elec_i,trial_i,dtw_i)  = templags(find(tempr == max(tempr)));
                        kdx = kdx+(windowlength/2);
                        dtw_i = dtw_i + 1;
                    else
                        brain_signal     = AllDataStructuresFT{sub_i,2}.trial{trial_i}(elec_i,kdx:end);
                        acoustic_signal  = envelope_music(trial_i,kdx:end);
                        [tempr,templags] = xcorr(zscore(brain_signal), ...
                                                 zscore(acoustic_signal),maxlag,'normalized');
                        rdata_music{sub_i}(elec_i,trial_i,dtw_i)    = max(tempr);
                        lagdata_music{sub_i}(elec_i,trial_i,dtw_i)  = templags(find(tempr == max(tempr)));
                        kdx = trialength;
                    end
                end
            end
        end
    else
        % croscorrelation if performed on the entire signal
        for trial_i=1:n_trials
            for elec_i=1:n_electrodes
                brain_signal     = AllDataStructuresFT{sub_i,2}.trial{trial_i}(elec_i,:);
                acoustic_signal  = envelope_music(trial_i,:);
                [tempr,templags] = xcorr(zscore(brain_signal), ...
                                         zscore(acoustic_signal),maxlag,'normalized');
                rdata_music{sub_i}(elec_i,trial_i,:)    = max(tempr);
                lagdata_music{sub_i}(elec_i,trial_i,:)  = templags(find(tempr == max(tempr)));
            end
        end
    end   
end

if segestimation == 1
    save([data_dir,filesep,'xcorr_' band2analyze '_windowed.mat'], ...
        'rdata_speech','lagdata_speech','rdata_music','lagdata_music', ...
        'band2analyze','sub2plot','AllChannelLabels','names4fields');
else
    save([data_dir,filesep,'xcorr_' band2analyze '.mat'], ...
        'rdata_speech','lagdata_speech','rdata_music','lagdata_music', ...
        'band2analyze','sub2plot','AllChannelLabels','names4fields');
end

% Permuted Croscorrelations
for sub_i=1:length(sub2plot)
    disp(['Estimating crosscorrelation for subject ' num2str(sub_i)]);
  
    for perm_i=1:nperms
        disp(['Permutation ' num2str(perm_i)]);
        
        RandTrialOrder = randperm(numel(AllDataStructuresFT{sub_i,1}.trial));
        
        if segestimation == 1
            for trial_i=1:RandTrialOrder
                for elec_i=1:size(AllDataStructuresFT{sub_i,1}.trial{trial_i},1)
                    
                    windowlength = 2 * fs; % use sliding windows of 2 seconds
                    kdx = 1; dtw_i = 1;                    
                    trialength = min([length(AllDataStructuresFT{sub_i,1}.trial{trial_i}(elec_i,:)) length(envelope_speech(trial_i,:))]);  %length of current trial

                    while kdx < trialength
                        if trialength - kdx > windowlength
                            [tempr,templags] = xcorr(zscore(AllDataStructuresFT{sub_i,1}.trial{trial_i}(elec_i,kdx:kdx+windowlength)),zscore(envelope_speech(RandTrialOrder(trial_i),kdx:kdx+windowlength)),maxlag,'normalized');
                            rdata_speech_perm{sub_i}(elec_i,trial_i,dtw_i,perm_i) = max(tempr);
                            lagdata_speech_perm{sub_i}(elec_i,trial_i,dtw_i,perm_i)= templags(find(tempr == max(tempr)));
                            [tempr,templags] = xcorr(zscore(AllDataStructuresFT{sub_i,2}.trial{trial_i}(elec_i,kdx:kdx+windowlength)),zscore(envelope_music(RandTrialOrder(trial_i),kdx:kdx+windowlength)),maxlag,'normalized');
                            rdata_music_perm{sub_i}(elec_i,trial_i,dtw_i,perm_i) = max(tempr);
                            lagdata_music_perm{sub_i}(elec_i,trial_i,dtw_i,perm_i)= templags(find(tempr == max(tempr)));
                            
                            kdx = kdx+(windowlength/2);
                            dtw_i = dtw_i + 1;
                        else
                            [tempr,templags] = xcorr(zscore(AllDataStructuresFT{sub_i,1}.trial{trial_i}(elec_i,kdx:end)),zscore(envelope_speech(RandTrialOrder(trial_i),kdx:end)),maxlag,'normalized');
                            rdata_speech_perm{sub_i}(elec_i,trial_i,dtw_i,perm_i) = max(tempr);
                            lagdata_speech_perm{sub_i}(elec_i,trial_i,dtw_i,perm_i)= templags(find(tempr == max(tempr)));
                            [tempr,templags] = xcorr(zscore(AllDataStructuresFT{sub_i,2}.trial{trial_i}(elec_i,kdx:end)),zscore(envelope_music(RandTrialOrder(trial_i),kdx:end)),maxlag,'normalized');
                            rdata_music_perm{sub_i}(elec_i,trial_i,dtw_i,perm_i) = max(tempr);
                            lagdata_music_perm{sub_i}(elec_i,trial_i,dtw_i,perm_i)= templags(find(tempr == max(tempr)));     
                            kdx = trialength;
                        end
                    end
                end
            end
        else            
            % PLV (using Florencia's function)
            for trial_i=1:RandTrialOrder
                for elec_i=1:size(AllDataStructuresFT{sub_i,1}.trial{trial_i},1)
                    [tempr,templags] = xcorr(zscore(AllDataStructuresFT{sub_i,1}.trial{trial_i}(elec_i,:)),zscore(envelope_speech(RandTrialOrder(trial_i),:)),maxlag,'normalized');
                    rdata_speech_perm{sub_i}(elec_i,trial_i,:,perm_i) = max(tempr);
                    lagdata_speech_perm{sub_i}(elec_i,trial_i,:,perm_i)= templags(find(tempr == max(tempr)));
                    [tempr,templags] = xcorr(zscore(AllDataStructuresFT{sub_i,2}.trial{trial_i}(elec_i,:)),zscore(envelope_music(RandTrialOrder(trial_i),:)),maxlag,'normalized');
                    rdata_music_perm{sub_i}(elec_i,trial_i,:,perm_i) = max(tempr);
                    lagdata_music_perm{sub_i}(elec_i,trial_i,:,perm_i)= templags(find(tempr == max(tempr)));
                end
            end
        end
    end
end

if segestimation == 1
     save(['E:\Matlab\IEEG\ECOGData_CROS_' band2analyze '_PERM2SEG.mat'],'rdata_speech_perm','lagdata_speech_perm','rdata_music_perm','lagdata_music_perm','band2analyze','sub2plot','AllChannelLabels','names4fields');
else
    save(['E:\Matlab\IEEG\ECOGData_CROS_' band2analyze '_PERM2.mat'],'rdata_speech_perm','lagdata_speech_perm','rdata_music_perm','lagdata_music_perm','band2analyze','sub2plot','AllChannelLabels','names4fields');
end
