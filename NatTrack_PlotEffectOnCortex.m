% NaturalisticTracking_ECOG project
%
% This script plots crosscorrelation coefficients and lags for the
% specified condition and frequency band. It also plots a histogram of the
% same values.

% S.Osorio - 2023

% initialize BST
% cd('E:\Matlab\brainstorm3');
% brainstorm
%%
clear, clc, close
condition2analyze = 'music'; % 'speech' or 'music'
band2analyze      = 'SFB';   % SFB (1-8 Hz) or HFB (70-120 Hz) 
segestimation     = 0;       % whether to use windowed data
plot_oncortex     = 1;       % whether to plot effect on cortical surface
plot_histogram    = 0;       % whether to plot histogram of plotted values
effect2plot       = 'lag';   % rho (corrcoefficients) or lag (xcorr lags)

% load data to plot
if segestimation == 1
    load(['E:\Matlab\IEEG\Data\CROSdata_',band2analyze,'_windowed.mat']);
else
    load(['E:\Matlab\IEEG\Data\CROSdata_',band2analyze,'.mat']);
end   

% set clustering parameters according to condition and frequency band of
% interest. These values were estimated in NatTrack_dbscan. Make sure they
% match the optimal parameters found in that script. 
if strcmpi(condition2analyze,'music') && strcmpi(band2analyze,'HFB')
        mindist = 0.018; minpoints = 12;
elseif strcmpi(condition2analyze,'speech') && strcmpi(band2analyze,'HFB')
        mindist = 0.010; minpoints = 12;
elseif strcmpi(condition2analyze,'both') && strcmpi(band2analyze,'HFB')
        mindist = 0.024; minpoints = 12;
elseif strcmpi(condition2analyze,'music') && strcmpi(band2analyze,'SFB')
        mindist = 0.018; minpoints = 12;
elseif strcmpi(condition2analyze,'speech') && strcmpi(band2analyze,'SFB')
        mindist = 0.018; minpoints = 12;
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
SurfaceFile   = 'E:\MATLAB\brainstorm_db\iEEG\anat\@default_subject\tess_cortex_pial_low.mat';

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
tmptable     = tabulate(testClusters)
deleteThis   = 0;

if any(tmptable(:,2) < minpoints)
    deleteThis = tmptable(find(tmptable(:,2) < minpoints),1);
end

% get rid of non-clustered electrodes
ValRange(testClusters == -1 | testClusters == deleteThis)     = [];
testMat(testClusters == -1 | testClusters == deleteThis,:,:)  = [];
subIDelec(testClusters == -1 | testClusters == deleteThis,:)  = [];
testClusters(testClusters == -1 | testClusters == deleteThis) = []; 
% plot effect on cortex
if plot_oncortex == 1
    [hFig, iDS, iFig] = view_surface(SurfaceFile);
    hFig.Color = [1 1 1];
    hold on;
    
    sh  = scatter3(testMat(:,1),testMat(:,2),testMat(:,3),60,ValRange,'filled');
    sh.MarkerFaceAlpha = .85;
    sh.SizeData = 100;

    ch = colorbar;
    if strcmpi(effect2plot,'rho')
        ch.Label.String = 'Rho';
        caxis([0.085 .13]);  % caxis values will likely need to be adjusted manually
    else
        ch.Label.String = 'Lag (ms)';
        caxis([-250 250]);
    end
    ch.FontSize = 20;
    title([condition2analyze ' - ' band2analyze],'FontWeight','normal','FontSize',22);
end

% histogram for ploted values
if plot_histogram == 1
    % colors per condition (1,:) music, (2,:) speech, (3,:) music and speech
    colors 	= [0.7176 0.2745 1.0000; ...
               0.9412 0.5804 0.3373; ...
               0.4667 0.6745 0.1882];  
           
    figure,clf
    hh = histogram(ValRange,5);
    if strcmpi(effect2plot,'rho')
        ylabel('Electrode count'); ylim([0 30]);
        xlabel('Correlation coefficient'); xlim([0 0.3]);
    else
        ylabel('Electrode count'); ylim([0 50]);
        xlabel('Lag (ms)'); xlim([-350 350]);
    end
    if strcmpi(condition2analyze,'music')
        hh.FaceColor = colors(1,:);
    else
        hh.FaceColor = colors(2,:);
    end       
    hh.EdgeColor = [1 1 1];
    title([condition2analyze ' - ' band2analyze],'FontWeight','Normal');
    set(gca,'FontSize',20,'FontName','Arial');
    box off
end