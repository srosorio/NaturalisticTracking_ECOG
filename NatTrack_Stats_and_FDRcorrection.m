% NaturalisticTracking_ECOG project
%
% This code estimates the critical correlation coefficient value corresponding  
% to a = 0.05 and performs FDR correction for observed p values for each
% subject. It also plots the total number of electrodes that survive
% statistics and saves a table with this same info at the subject level.

% S.Osorio - 2023

clear, clc,

% create a directory to save data
iEEG_dir = 'E:\Matlab\IEEG';
data_dir = [iEEG_dir,filesep,'Data'];

% SFB (1-8 Hz) or HFB (70-120 Hz)
band2analyze  = 'HFB';
% use sliding window data (1 = yes, 0 = no)
segestimation = 0;
% plot histogram of FDR-corrected p values (1 = yes, 0 = no)
plot_FDR      = 0;
% plot number of electrodes per condition (1 = yes, 0 = no)
plot_bars     = 1;

% load data
if segestimation == 1
    load([data_dir,filesep,'xcorr_',band2analyze,'_windowed.mat'])
    load([data_dir,filesep,'xcorr_',band2analyze,'_windowed_PERM.mat'])
else
    load([data_dir,filesep,'xcorr_',band2analyze,'.mat'])
    load([data_dir,filesep,'xcorr_',band2analyze,'_PERM.mat'])
end

alpha       = 0.05;  % p value       
FDRcorrect  = 0;     % do FDR correction?   
n_trials    = size(r_speech{1},2);
n_conds     = 2;

% this is to merge High Density data to their corresponding participants
for sub_i=1:length(sub2plot)
    if contains(sub2plot{sub_i},'HD')
        r_music{sub_i-1}       = [r_music{sub_i-1}; r_music{sub_i}];
        r_music_perm{sub_i-1}  = [r_music_perm{sub_i-1}; r_music_perm{sub_i}];
        lag_music{sub_i-1}     = [lag_music{sub_i-1}; lag_music{sub_i}];
        r_speech{sub_i-1}      = [r_speech{sub_i-1}; r_speech{sub_i}];
        r_speech_perm{sub_i-1} = [r_speech_perm{sub_i-1}; r_speech_perm{sub_i}];
        lag_speech{sub_i-1}    = [lag_speech{sub_i-1}; lag_speech{sub_i}];
        
        AllChannelLabels{sub_i-1} = [AllChannelLabels{sub_i-1};  AllChannelLabels{sub_i}];
    end
end
% now delete redundant data
r_music(contains(sub2plot,'HD'))        = [];
r_music_perm(contains(sub2plot,'HD'))   = [];
lag_music(contains(sub2plot,'HD'))      = [];
r_speech(contains(sub2plot,'HD'))       = [];
r_speech_perm(contains(sub2plot,'HD'))  = [];
lag_speech(contains(sub2plot,'HD'))     = [];
AllChannelLabels(contains(sub2plot,'HD')) = [];
sub2plot(contains(sub2plot,'HD'))       = [];

% initialize some data containers we will need for later
[pvals_speech, pvals_music, lags_speech, lags_music, fdr_speech, fdr_music] = deal(cell(1,length(sub2plot)));
n_subs  = size(sub2plot,2);

% loop through subjects and electrodes to get corresponding p values
for sub_i=1:n_subs
    n_elecs  = size(r_speech{sub_i},1);
    n_perms  = size(r_speech_perm{sub_i},4);
    % obtain the mean p value across trials for each electrode within each subject
    for trial_i=1:n_trials           
        for elec_i=1:n_elecs
            % for analysis conducted in moving windows
            if segestimation == 1   
                % raw p values for speech
                % null distribution from permuted data
                null_dist   = mean(mean(r_speech_perm{sub_i}(elec_i,:,:,:),3),2);
                % observed correlation coefficient
                observed_r  = mean(mean(r_speech{sub_i}(elec_i,:,:),3));
                pvals_speech{sub_i}(elec_i,:) = sum(null_dist>=observed_r) /  n_perms;             
                % raw p values for music
                null_dist   = mean(mean(r_music_perm{sub_i}(elec_i,:,:,:),3),2);
                observed_r  = mean(mean(r_music{sub_i}(elec_i,:,:),3));   
                pvals_music{sub_i}(elec_i,:)  = sum(null_dist>=observed_r) /  n_perms;
            % for analysis conducted in the entire signal
            else
                % speech
                null_dist   = mean(squeeze(r_speech_perm{sub_i}(elec_i,:,:,:)),1);
                observed_r  = mean(r_speech{sub_i}(elec_i,:,:),2);
                pvals_speech{sub_i}(elec_i,:) = sum(null_dist >= observed_r) / n_perms;
                % music 
                null_dist   = mean(squeeze(r_music_perm{sub_i}(elec_i,:,:,:)),1);
                observed_r  = mean(r_music{sub_i}(elec_i,:,:),2);
                pvals_music{sub_i}(elec_i,:)  = sum(null_dist >= observed_r) / n_perms;
            end
        end
    end
    % p val histograms
    if plot_FDR == 1
        figure(1)
        subplot(6,5,sub_i)
        histogram(pvals_speech{sub_i})
        xlim([0 1]); ylim([0 30])
        figure(2)
        subplot(6,5,sub_i)
        xlim([0 1]); ylim([0 30])
        histogram(pvals_music{sub_i})
    end
	% implement FDR correction (this is done within each subject)
    if FDRcorrect == 1
        fdr_speech{sub_i} = mafdr(pvals_speech{sub_i},'lambda',0.15);
        fdr_music{sub_i}  = mafdr(pvals_music{sub_i},'lambda',0.15);
    else
        fdr_speech{sub_i} = pvals_speech{sub_i};
        fdr_music{sub_i}  = pvals_music{sub_i};
    end
    % FDR-corrected p val histograms
    if plot_FDR == 1
        figure(3)
        subplot(6,5,sub_i)
        histogram(fdr_speech{sub_i})
        xlim([0 1]); ylim([0 30])    
        figure(4)
        subplot(6,5,sub_i)
        histogram(fdr_music{sub_i})
        xlim([0 1]); ylim([0 30])  
    end

    % finally, obtain the mean lag (ms) values across trials per electrode
    lags_speech{sub_i} = mean(mean(lag_speech{sub_i},3),2);
    lags_music{sub_i}  = mean(mean(lag_music{sub_i},3),2);
end

for sub_i=1:n_subs
    % find number of electrodes per subject
    matlength(sub_i) = numel(fdr_speech{sub_i});    
end

TotalElecs = sum(matlength); % total electrodes
maxnelecs  = max(matlength); % max number of electrodes across subjects

% now create a numeric array for all p values and another one for all lags 
PValMat    = NaN(maxnelecs,n_subs,n_conds); 
LagMat     = NaN(maxnelecs,n_subs,n_conds);

% store data for all subjects 
for sub_i=1:n_subs
    clear tmpSpeech tmpMusic
    tmpSpeech = fdr_speech{sub_i};
    tmpMusic  = fdr_music{sub_i};

    PValMat(1:length(tmpSpeech),sub_i,1) = tmpSpeech;
    LagMat(1:length(tmpSpeech),sub_i,1)  = lags_speech{sub_i};
    PValMat(1:length(tmpMusic),sub_i,2)  = tmpMusic;    
    LagMat(1:length(tmpSpeech),sub_i,2)  = lags_music{sub_i};    
end
 
% now, we look for the correlation coefficients where p < alpha and store them
for sub_i=1:n_subs
    r_data_speech = NaN(maxnelecs,n_trials);
    r_data_music  = NaN(maxnelecs,n_trials);   
    for trial_i=1:n_trials
        n_elecs = size(r_speech{sub_i},1);
        for p_i=1:n_elecs
%             disp(['Calculating ?Coh for electrode ' num2str(p_i) ' for subject ' num2str(sub_i)]);
            if fdr_speech{sub_i}(p_i) < alpha
                r_data_speech(p_i,trial_i) = mean(r_speech{sub_i}(p_i,trial_i,:),3);
            end
            if fdr_music{sub_i}(p_i) < alpha
                r_data_music(p_i,trial_i) = mean(r_music{sub_i}(p_i,trial_i,:),3);
            end
        end
    end    
    % and now we get the mean correlation coefficient across trials
    r_data{sub_i,1} = nanmean(r_data_speech,2);
    r_data{sub_i,2} = nanmean(r_data_music,2);
end

% finally, let's create an array to store all statistically significant
% correlation values (:,:,1) and lag values (:,:,2) for all subjects (dim
% 2) and for each subjects electrodes (dim 1)
dataMat   = NaN(maxnelecs,n_subs,n_conds);

for sub_i=1:size(r_data,1)
    clear tmpSpeech tmpMusic
    tmpSpeech = r_data{sub_i,1};
    tmpMusic  = r_data{sub_i,2};

    dataMat(1:length(tmpSpeech),sub_i,1) = tmpSpeech;
    LagMat(isnan(r_data{sub_i,1}),sub_i,1) = NaN;
    PValMat(isnan(r_data{sub_i,1}),sub_i,1) = NaN;
    dataMat(1:length(tmpMusic),sub_i,2)  = tmpMusic;    
    LagMat(isnan(r_data{sub_i,2}),sub_i,2) = NaN;
    PValMat(isnan(r_data{sub_i,2}),sub_i,2) = NaN;
end

% print how many subjects show a statistically significant effect per condition
subeffect   = sum(any(dataMat(:,:,1)));
disp(['SPEECH: Statistically significant data in ' num2str(subeffect) ' out of ' num2str(length(sub2plot)) ' subjects'])
subeffect   = sum(any(dataMat(:,:,2)));
disp(['MUSIC: Statistically significant data in ' num2str(subeffect) ' out of ' num2str(length(sub2plot)) ' subjects'])

% save data
if segestimation == 1
    save([data_dir,filesep,'CROSdata_' band2analyze '_windowed.mat'],  ... 
    'dataMat','LagMat','PValMat','sub2plot','TotalElecs','AllChannelLabels','names4fields');
else
    save([data_dir,filesep,'CROSdata_' band2analyze '.mat'],  ...
        'dataMat','LagMat','PValMat','sub2plot','TotalElecs','AllChannelLabels','names4fields');
end

% bar plots for total number of electrodes per condition
if plot_bars == 1

    ElecsperSub_Speech = sum(~isnan(dataMat(:,:,1)));
    ElecsperSub_Music  = sum(~isnan(dataMat(:,:,2)));
    ElecsperSub_Both   = sum(~isnan(dataMat(:,:,1)) & ~isnan(dataMat(:,:,2)));
    
    % save data
    save([data_dir,filesep,'ElecsPerSub_',band2analyze],'ElecsperSub_Speech','ElecsperSub_Music','ElecsperSub_Both')
    
    % create bar plot
    colors = [0.7176 0.2745 1.0000; ...
              1.0000 0.4118 0.1608; ...
              0.4667 0.6745 0.1882];    
    ElecsPerSub = [ElecsperSub_Music',ElecsperSub_Speech',ElecsperSub_Both'];
    
    % let's create a table with basic statistics
    RhosMusic  = dataMat(:,:,2);
    RhosMusic  = RhosMusic(~isnan(RhosMusic));
    RhosSpeech = dataMat(:,:,1);
    RhosSpeech = RhosSpeech(~isnan(RhosSpeech));
    pMusic     = PValMat(:,:,2);
    pMusic     = pMusic(~isnan(dataMat(:,:,2)));
    pSpeech    = PValMat(:,:,1);
    pSpeech    = pSpeech(~isnan(dataMat(:,:,1)));
    RhosBoth   = dataMat(~isnan(dataMat));
    pBoth      = PValMat(~isnan(dataMat));
    
    disp(['N music = ' num2str(length(RhosMusic))]);
    disp(['N speech = ' num2str(length(RhosSpeech))]);

    StatsTable = table([mean(RhosMusic); std(RhosMusic); mean(pMusic)], ...
                       [mean(RhosSpeech);std(RhosSpeech);mean(pSpeech)], ...
                       [mean(RhosBoth); std(RhosBoth); mean(pBoth)], ...
                       'VariableNames',{'Music','Speech','Both'}, ...
                       'RowNames',{'mean rho','SD rho','mean p'});
    display(StatsTable);
    
    figure(1), clf
    bh1 = bar(sum(ElecsPerSub),'FaceColor','flat');
    bh1.CData = colors;
    bh1.EdgeColor = [1 1 1];
    ylim([0 250]);
    xticklabels({'Music','Speech','Both'});
    xtickangle(45);
    ylabel('Electrode count');
    set(gca,'FontSize',20,'FontName','Arial');
    title(band2analyze,'FontWeight','normal');
    box off
    
    % save number of electrodes per subject as a table
    ElecsPersSub = array2table(ElecsPerSub,'RowNames',sub2plot,'VariableNames',{'Music','Speech','Both'});
    writetable(ElecsPersSub,[data_dir,filesep,'ElecsPerSub_',band2analyze,'.xlsx'])
end