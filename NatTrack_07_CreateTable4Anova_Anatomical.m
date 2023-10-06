clear, clc

% set paths
iEEG_dir = 'F:\Matlab\IEEG';
data_dir = [iEEG_dir,filesep,'Data'];
cd(data_dir)

fnames          = dir('Electrode_AnatLabels*');
print_diagnosis = 1;

if print_diagnosis == 1
    for idx=1:length(fnames)
        load(fnames(idx).name)
        AnatLabels = cellstr(AnatLabels);
        tmp1 = contains(AnatLabels,'premotor') | contains(AnatLabels,'somatomotor') | contains(AnatLabels,'somatosensory');
        %     tmp2 = contains(AnatLabels,'stg') | contains(AnatLabels,'mtg');
        AnatLabels(tmp1) = {'sensorymotor'};
        %     AnatLabels(tmp2) = {'temporal'};
        info = split(fnames(idx).name,'_');
        cond = info{3};
        band = split(info{4},'.mat');
        band = band{1};
        disp([cond, '-' ,band]);
        tbl = sortrows(tabulate(cellstr(AnatLabels)),3,'descend');
        disp(tbl);
    end
end

condition2analyze = {'Speech','Music'};   % 'speech' or 'music'
band2analyze      = {'SFB','HFB'};      % SFB (1-8 Hz) or HFB (70-120 Hz) 
effect2plot       = {'rho','lag'};      % rho (corrcoefficients) or lag (xcorr lags)
perm_type         = 'wn';

% set paths
iEEG_dir = 'F:\Matlab\IEEG';
data_dir = [iEEG_dir,filesep,'Data'];

for band_i=1:length(band2analyze)
    for cond_i=1:length(condition2analyze)
        load(['F:\Matlab\IEEG\Data\CROSdata_',band2analyze{band_i},'_whitenoise.mat']);
        
        %if using whitenoise permutations
        if strcmpi(condition2analyze{cond_i},'speech') & strcmpi(band2analyze{band_i},'SFB')
            mindist = 0.006; minpoints = 16;
        elseif strcmpi(condition2analyze{cond_i},'speech') & strcmpi(band2analyze{band_i},'HFB')
            mindist = 0.006; minpoints = 12;
        elseif strcmpi(condition2analyze{cond_i},'music') & strcmpi(band2analyze{band_i},'SFB')
            mindist = 0.012; minpoints = 12;
        elseif strcmpi(condition2analyze{cond_i},'music') & strcmpi(band2analyze{band_i},'HFB')
            mindist = 0.014; minpoints = 12;
        end
        
        % locate data in separate variables and plot
        rho4speech = dataMat(:,:,1);
        rho4music  = dataMat(:,:,2);
        lag4speech = LagMat(:,:,1);
        lag4music  = LagMat(:,:,2);
        
        % create arrays containing the data of interest
        if strcmpi(condition2analyze{cond_i},'speech')
            counter = 1;
            for sub_i=1:length(sub2plot)
                clear ThisSubStruct
                ThisSubStruct = cell2struct(AllChannelLabels{sub_i},names4fields,2);
                for idx=1:length(dataMat)
                    if ~isnan(rho4speech(idx,sub_i)) %% && isnan(rhos4music(idx,sub_i)) %%
                        testMat(counter,:)   = [ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3)];
                        RhoRange(counter,:)  = rho4speech(idx,sub_i);
                        LagRange(counter,:)  = lag4speech(idx,sub_i);
                        subIDelec(counter,:) = [sub_i,idx];
                        counter = counter + 1;
                    end
                end
            end
        elseif strcmpi(condition2analyze{cond_i},'music')
            counter = 1;
            for sub_i=1:length(sub2plot)
                clear ThisSubStruct
                ThisSubStruct = cell2struct(AllChannelLabels{sub_i},names4fields,2);
                for idx=1:length(dataMat)
                    if  ~isnan(rho4music(idx,sub_i)) %% isnan(rhos4speech(idx,sub_i)) &&
                        testMat(counter,:)   = [ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3)];
                        RhoRange(counter,:)  = rho4music(idx,sub_i);
                        LagRange(counter,:)  = lag4music(idx,sub_i);
                        subIDelec(counter,:) = [sub_i,idx];
                        counter = counter + 1;
                    end
                end
            end
        end
        
        % now get only the desired electrodes that belong to the specified cluster
        testClusters = dbscan(testMat,mindist,minpoints);
        tmptable     = sortrows(tabulate(testClusters),-3)
        deleteThis   = 0;
        
        % identify clusters where number of elecs < minpoints
        if any(tmptable(:,2) < minpoints)
            deleteThis = tmptable(tmptable(:,2) < minpoints,1);
        end
        
        % group smaller surviving clusters into one single cluster
        if size(tmptable,1)-1 > 3
            testClusters(ismember(testClusters,tmptable(4:end,1))) = 3;
        end
        
        % get rid of non-clustered electrodes
        RhoRange(ismember(testClusters,[-1; deleteThis]))     = [];
        LagRange(ismember(testClusters,[-1; deleteThis]))     = [];
        testMat(ismember(testClusters,[-1; deleteThis]),:,:)  = [];
        subIDelec(ismember(testClusters,[-1; deleteThis]),:)  = [];
%         testClusters(ismember(testClusters,[-1; deleteThis])) = [];

        % get anatomical labels to create anova table
        load([data_dir,filesep,'Electrode_AnatLabels_' condition2analyze{cond_i} '_' band2analyze{band_i} '.mat'])
        AnatLabels = cellstr(AnatLabels);
        tmp1 = contains(AnatLabels,'premotor') | contains(AnatLabels,'somatomotor') | contains(AnatLabels,'somatosensory');
        %     tmp2 = contains(AnatLabels,'stg') | contains(AnatLabels,'mtg');
        AnatLabels(tmp1) = {'sensorymotor'};
        %     AnatLabels(tmp2) = {'temporal'};
        disp([condition2analyze{cond_i}, '-' ,band2analyze{band_i}]);
        tbl = sortrows(tabulate(cellstr(AnatLabels)),3,'descend');
        disp(tbl);
        % get only the data we need
        these_elecs = contains(AnatLabels,'mtg') | ...
            contains(AnatLabels,'stg') | ...
            contains(AnatLabels,'supramarginal') | ...
            contains(AnatLabels,'sensorymotor');
        % create and save table for this effect
        if band_i==1 && cond_i == 1
            tbl4anova = table(subIDelec(these_elecs,1), ...
                    repmat(convertCharsToStrings(band2analyze{band_i}),sum(these_elecs),1), ...
                    repmat(convertCharsToStrings(condition2analyze{cond_i}),sum(these_elecs),1), ...
                    AnatLabels(these_elecs)', ...
                    RhoRange(these_elecs), LagRange(these_elecs), 'VariableNames', {'subject','frequency','condition','region','r','lag'});
        else
             newtbl4anova = table(subIDelec(these_elecs,1), ...
                    repmat(convertCharsToStrings(band2analyze{band_i}),sum(these_elecs),1), ...
                    repmat(convertCharsToStrings(condition2analyze{cond_i}),sum(these_elecs),1), ...
                    AnatLabels(these_elecs)', ...
                    RhoRange(these_elecs), LagRange(these_elecs), 'VariableNames', {'subject','frequency','condition','region','r','lag'});           
             tbl4anova = [tbl4anova; newtbl4anova];
        end
    end
end

% save data
writetable(tbl4anova,[data_dir,filesep,'AnovaNatTrackTable_wn_ANAT.csv'])
