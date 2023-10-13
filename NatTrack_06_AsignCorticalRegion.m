% NaturalisticTracking_ECOG project
%
% This script plots crosscorrelation coefficients and lags for the
% specified condition and frequency band. It also plots a histogram of the
% same values.

% S.Osorio - 2023

% initialize BST
% cd('F:\Matlab\brainstorm3');
% brainstorm
%%
clear, clc, close all
condition2analyze = 'Music';   % 'speech' or 'music'
band2analyze      = 'HFB';      % SFB (1-8 Hz) or HFB (70-120 Hz) 
segestimation     = 0;          % whether to use windowed data
plot_oncortex     = 1;          % whether to plot effect on cortical surface
plot_histogram    = 1;          % whether to plot histogram of plotted values
effect2plot       = 'lag';      % rho (corrcoefficients) or lag (xcorr lags)
perm_type         = 'wn';

% set paths
iEEG_dir = 'F:\Matlab\IEEG';
data_dir = [iEEG_dir,filesep,'Data'];

% load data to plot
if segestimation == 1
    load(['F:\Matlab\IEEG\Data\CROSdata_',band2analyze,'_windowed.mat']);
elseif strcmpi(perm_type,'ts')
    load(['F:\Matlab\IEEG\Data\CROSdata_',band2analyze,'_trialshuffle.mat']);
elseif strcmpi(perm_type,'wn')
    load(['F:\Matlab\IEEG\Data\CROSdata_',band2analyze,'_whitenoise.mat']);
end   

% these are the parameters after optimization per condition and freq band
if strcmpi(perm_type,'ts') 
    % if using trial shuffling permutations
    if strcmpi(condition2analyze,'speech') && strcmpi(band2analyze,'SFB')
        mindist = 0.016; minpoints = 14;
    elseif strcmpi(condition2analyze,'speech') && strcmpi(band2analyze,'HFB')
        mindist = 0.01; minpoints = 12;
    elseif strcmpi(condition2analyze,'music') && strcmpi(band2analyze,'SFB')
        mindist = 0.02; minpoints = 14;
    elseif strcmpi(condition2analyze,'music') && strcmpi(band2analyze,'HFB')
        mindist = 0.018; minpoints = 12;        
    end
else
    %if using whitenoise permutations
    if strcmpi(condition2analyze,'speech') && strcmpi(band2analyze,'SFB')
        mindist = 0.006; minpoints = 16;
    elseif strcmpi(condition2analyze,'speech') && strcmpi(band2analyze,'HFB')
        mindist = 0.006; minpoints = 12;
    elseif strcmpi(condition2analyze,'music') && strcmpi(band2analyze,'SFB')
        mindist = 0.012; minpoints = 12;
    elseif strcmpi(condition2analyze,'music') && strcmpi(band2analyze,'HFB')
        mindist = 0.014; minpoints = 12;
    end
end

% locate data in separate variables and plot
if strcmpi(effect2plot,'rho')
    effect4speech = dataMat(:,:,1);
    effect4music  = dataMat(:,:,2);
elseif strcmpi(effect2plot,'lag')
    effect4speech = LagMat(:,:,1);
    effect4music  = LagMat(:,:,2);
end

% get MNI cortical surface using brainstorm. We will plot data here.
SurfaceFile   = 'F:\MATLAB\brainstorm_db\iEEG\anat\@default_subject\tess_cortex_pial_low.mat';

% create arrays containing the data of interest
if strcmpi(condition2analyze,'speech')
    counter = 1;
    for sub_i=1:length(sub2plot)    
        clear ThisSubStruct
        ThisSubStruct = cell2struct(AllChannelLabels{sub_i},names4fields,2);        
        for idx=1:length(dataMat)
            if ~isnan(effect4speech(idx,sub_i)) %% && isnan(rhos4music(idx,sub_i)) %%
                testMat(counter,:)   = [ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3)];
                ValRange(counter,:)  = effect4speech(idx,sub_i);
                subIDelec(counter,:) = [sub_i,idx];
                counter = counter + 1;
            end
        end
    end    
elseif strcmpi(condition2analyze,'music')    
    counter = 1;    
    for sub_i=1:length(sub2plot)       
        clear ThisSubStruct
        ThisSubStruct = cell2struct(AllChannelLabels{sub_i},names4fields,2);        
        for idx=1:length(dataMat)
            if  ~isnan(effect4music(idx,sub_i)) %% isnan(rhos4speech(idx,sub_i)) &&
                testMat(counter,:)   = [ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3)];
                ValRange(counter,:)  = effect4music(idx,sub_i);
                subIDelec(counter,:) = [sub_i,idx];
                counter = counter + 1;
            end
        end
    end    
end

% now get only the desired electrodes that belong to the specified cluster
testClusters = dbscan(testMat,mindist,minpoints);
tmptable = sortrows(tabulate(testClusters),-3)
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
ValRange(ismember(testClusters,[-1; deleteThis]))     = [];
testMat(ismember(testClusters,[-1; deleteThis]),:,:)  = [];
subIDelec(ismember(testClusters,[-1; deleteThis]),:)  = [];
testClusters(ismember(testClusters,[-1; deleteThis])) = []; 

% get mean rho vals at the subject level
if strcmpi(condition2analyze,'music')
    for i=1:30
        meanRho(i)  = mean(dataMat(subIDelec(subIDelec(:,1) == i,2),i,2));
        meanPval(i) = mean(PValMat(subIDelec(subIDelec(:,1) == i,2),i,2));
        meanLag(i)  = mean(LagMat(subIDelec(subIDelec(:,1) == i,2),i,2));
    end
else 
    for i=1:30
        meanRho(i)  = mean(dataMat(subIDelec(subIDelec(:,1) == i,2),i,1));
        meanPval(i) = mean(PValMat(subIDelec(subIDelec(:,1) == i,2),i,1));
        meanLag(i)  = mean(LagMat(subIDelec(subIDelec(:,1) == i,2),i,1));        
    end
end

% save data
% save([data_dir,filesep,'stats_clustered_',band2analyze,'_',condition2analyze,'.mat'],'meanRho','meanPval','meanLag');

% plot statistical effect on cortex
if plot_oncortex == 1

    for jdx=1:length(testMat)
        close
        [hFig, iDS, iFig] = view_surface(SurfaceFile);
        hFig.Color = [1 1 1];
        hold on;
        sh = scatter3(testMat(jdx,1),testMat(jdx,2),testMat(jdx,3),60,'filled');
        tmp = input(['Where is the electrode ' num2str(jdx) '? >>> '],'s');
        AnatLabels{jdx} = convertCharsToStrings(tmp);
    end
end

 save([data_dir,filesep,'Electrode_AnatLabels_' condition2analyze, '_', band2analyze '.mat'],'AnatLabels')
%     sh  = scatter3(testMat(:,1),testMat(:,2),testMat(:,3),60,ValRange,'filled');
%     sh.MarkerFaceAlpha = .85;
%     sh.SizeData = 100;

%     ch = colorbar;
%     if strcmpi(effect2plot,'rho')
%         ch.Label.String = 'Rho';
%         %caxis([0.05 .10]);  % caxis values will likely need to be adjusted manually
%     else
%         ch.Label.String = 'Lag (ms)';
%         caxis([-250 250]);
%     end
%     ch.FontSize = 20;
%     title([condition2analyze ' - ' band2analyze],'FontWeight','normal','FontSize',22);
% end
% 
% % histogram for ploted values
% if plot_histogram == 1
%     % colors per condition (1,:) music, (2,:) speech, (3,:) music and speech
%     colors 	= [0.7176 0.2745 1.0000; ...
%                0.9412 0.5804 0.3373; ...
%                0.4667 0.6745 0.1882];  
%            
%     figure,clf
%     hh = histogram(ValRange);
%     if strcmpi(effect2plot,'rho')
%         ylabel('Electrode count'); %ylim([0 30]);
%         xlabel('Correlation coefficient'); xlim([0 0.3]);
%     else
%         ylabel('Electrode count'); %ylim([0 50]);
%         xlabel('Lag (ms)'); xlim([-350 350]);
%     end
%     if strcmpi(condition2analyze,'music')
%         hh.FaceColor = colors(1,:);
%     else
%         hh.FaceColor = colors(2,:);
%     end       
%     hh.EdgeColor = [1 1 1];
%     title([condition2analyze ' - ' band2analyze],'FontWeight','Normal');
%     set(gca,'FontSize',20,'FontName','Arial');
%     box off
% end
% 
% % print some statistics
% if strcmpi(condition2analyze,'music')
%     pvals = PValMat(:,:,2);
% else
%     pvals = PValMat(:,:,1);
% end
% pvals = pvals(~isnan(pvals));
% 
% if strcmpi(effect2plot,'rho')
%     disp(['n = ' num2str(sum(~isnan(meanRho))), ', r = ' num2str(mean(ValRange),2), ', elecs = ' num2str(length(ValRange)), ', p = ' num2str(mean(pvals),2), ', SD = ' num2str(std(ValRange),2)]);
% else
%     disp(['n = ' num2str(sum(~isnan(meanRho))), ', r = ' num2str(mean(ValRange),2), ', elecs = ' num2str(length(ValRange)), ', p = ' num2str(mean(pvals),2), ', SD = ' num2str(std(ValRange),2)]);
% end