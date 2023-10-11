% NaturalisticTracking_ECOG project
%
% This script plots statistically significant electrodes in common space
% (MNI) per condition, and conducts a data-driven cluster analysis to identify
% cortical regions driving statistical effects.

% S.Osorio - 2023

% initialize brainstorm (for data visualization)
cd('F:\Matlab\brainstorm3');
brainstorm
%%
clear, clc, close

% -------------------------------------------------------------------------
%                             SET PARAMETERS
% -------------------------------------------------------------------------
segestimation     = 0;  % wether to use segmented data
plot_neteffect    = 0;  % plot effect prior to statistics
dbscan_param      = 0;  % plot dbscan optimization process
plot_clusters     = 1;  % only after dbscan optimization
condition2analyze = 'both';
band2analyze      = 'SFB';
perm_type         = 'wn';   % wn = whitenoise, ts = trialshuffling
% -------------------------------------------------------------------------

% set paths
iEEG_dir = 'F:\Matlab\IEEG';
data_dir = [iEEG_dir,filesep,'Data'];

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

% ------------------- nothing to change from here on ----------------------
% mindist = 0.004; minpoints = 12;

% colors per condition (1,:) music, (2,:) speech, (3,:) music and speech
colors  = [0.7176 0.2745 1.0000; ...
           0.9412 0.5804 0.3373; ...
           0.4667 0.6745 0.1882];              
    
% load data to plot
if segestimation == 1
    load(['F:\Matlab\IEEG\Data\CROSdata_',band2analyze,'_windowed.mat']);
elseif strcmpi(perm_type,'ts')
    load(['F:\Matlab\IEEG\Data\CROSdata_',band2analyze,'_trialshuffle.mat']);
elseif strcmpi(perm_type,'wn')
    load(['F:\Matlab\IEEG\Data\CROSdata_',band2analyze,'_whitenoise.mat']);
end    

% locate data in separate variables and plot
rhos4speech = dataMat(:,:,1);
rhos4music  = dataMat(:,:,2);

% get MNI cortical surface using brainstorm. We will plot data here.
SurfaceFile   = 'F:\MATLAB\brainstorm_db\iEEG\anat\@default_subject\tess_cortex_pial_low.mat'; %['C:\Users\andre\OneDrive\Documentos\MATLAB\brainstorm_db\iEEG\anat\' sub2plot '\tess_cortex_central_low.mat'];

n_subs  = length(sub2plot);
n_elecs = length(dataMat);

% first, we get the info we need and put it in a single array 
if strcmpi(condition2analyze,'speech')
    counter = 1;
    for sub_i=1:n_subs    
        % get the subject structure with electrode info
        clear ThisSubStruct
        ThisSubStruct = cell2struct(AllChannelLabels{sub_i},names4fields,2);    
        for idx=1:n_elecs
            % electrodes selective to speech 
            if ~isnan(rhos4speech(idx,sub_i))
                datamat(counter,:)   = [ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3)];
                ValRange(counter,:)  = rhos4speech(idx,sub_i);
                subIDelec(counter,:) = [sub_i,idx];
                counter = counter + 1;
            end
        end
    end    
elseif strcmpi(condition2analyze,'music')    
    counter = 1;   
    for sub_i=1:n_subs       
        clear ThisSubStruct
        ThisSubStruct = cell2struct(AllChannelLabels{sub_i},names4fields,2);
        for idx=1:n_elecs
            % electrodes selective to music
            if  ~isnan(rhos4music(idx,sub_i))
                datamat(counter,:)   = [ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3)];
                ValRange(counter,:)  = rhos4speech(idx,sub_i);
                subIDelec(counter,:) = [sub_i,idx];
                counter = counter + 1;
            end
        end
    end 
elseif strcmpi(condition2analyze,'both')
    counter = 1;   
    music_info = load([data_dir,filesep,'dbscan_results_music_' band2analyze '.mat']);
    speech_info = load([data_dir,filesep,'dbscan_results_speech_' band2analyze '.mat']);
    for sub_i=1:n_subs  
        tmp = [music_info.subIDelec(music_info.subIDelec(:,1) == sub_i,2); speech_info.subIDelec(speech_info.subIDelec(:,1) == sub_i,2)];
        [gc,grps] = groupcounts(tmp); clear tmp;
        shared_elecs = grps(gc == 2);
        clear ThisSubStruct
        ThisSubStruct = cell2struct(AllChannelLabels{sub_i}(shared_elecs,:),names4fields,2);    
        for idx=1:length({ThisSubStruct.Name})
            % electrodes selective to bot speech and music
%             if  ~isnan(rhos4music(idx,sub_i)) && ~isnan(rhos4speech(idx,sub_i)) 
                datamat(counter,:)   = [ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3)];
                ValRange(counter,:)  = rhos4speech(idx,sub_i);
                subIDelec(counter,:) = [sub_i,idx];
                counter = counter + 1;
%             end
        end
    end  
end

% plot all electrodes (Net effect, without cluster analysis)
if plot_neteffect == 1
    % bst will throw an error every time the cortical surface is ploted for
    % the first time. This try catch statement overrides this error. Close
    % the loading bar if it doesn't disappear on its own. 
    try
        [hFig, iDS, iFig] = view_surface(SurfaceFile);
    catch
        [hFig, iDS, iFig] = view_surface(SurfaceFile);
    end
    hFig.Color = [1 1 1];
    hold on;
    % plot electrodes
    sh  = scatter3(datamat(:,1),datamat(:,2),datamat(:,3),60,'filled');
    sh.MarkerFaceAlpha = .8;
    sh.SizeData = 100;
    % set colors according to condition
    if strcmpi(condition2analyze,'music')
        sh.CData = colors(1,:);
    elseif strcmpi(condition2analyze,'speech')
        sh.CData = colors(2,:);
    elseif strcmpi(condition2analyze,'both')
        sh.CData = colors(3,:);
    end
end

% this is to determine the best parameters for dbscan algorithm. We are
% looking for the highest possible value in the y axis and the lowest
% possible value in the x axis
if dbscan_param == 1
    clear testClusters prctelecs numclusters numoutliers
    dsts   = 0:0.002:0.04;
    points = 12:2:20;

    for idx=1:length(dsts)
        clear clusterdata
        for jdx=1:length(points)
            tmp = dbscan(datamat,dsts(idx),points(jdx));
            clusterdata = tabulate(tmp);
            if numel(clusterdata(:,1)) == 1 && (clusterdata(1,1) == -1)
                prctelecs(jdx,idx)   = 0;
                prctoutliers(jdx,idx) = clusterdata(1,3);
            elseif numel(clusterdata(:,1)) > 1 && (clusterdata(1,1) == -1)
                prctelecs(jdx,idx)   = sum(clusterdata(2:end,3));
                numclusters(jdx,idx) = sum(clusterdata(:,1) ~= -1);
                prctoutliers(jdx,idx) = clusterdata(1,3);
            elseif numel(clusterdata(:,1)) == 1 && (clusterdata(1,1) == 1)
                prctelecs(jdx,idx)   = sum(clusterdata(:,3));
                numclusters(jdx,idx) = sum(clusterdata(:,1) ~= -1);
                prctoutliers(jdx,idx) = 0;
            end
        end
    end
    
    opt_idx = prctelecs ./ (dsts*1000) .* points' .* numclusters;  
    max_val = max(max(opt_idx));
    opt_idx = opt_idx / max_val; 
    
    % plot data
    figure(3), clf
    plot(dsts,opt_idx,'LineWidth',1.5)
    ylabel('OI_{norm}');
    xlabel('epsilon');
    legend(num2str(points'),'Location','SouthEast'); legend boxoff
    ylim([0 1.02]); set(gca,'FontSize',20);
    title([condition2analyze ' - ' band2analyze],'FontWeight','normal')
    box off
end

% finally, this is to plot the clusters using the best parameters as
% determined above
if plot_clusters == 1
    %colors for different clusters within conditions
    if strcmpi(condition2analyze,'music')
        colors = [0.7176 0.2745 1.0000];
    elseif strcmpi(condition2analyze,'speech')
        colors = [1.0000 0.4118 0.1608];
    elseif strcmpi(condition2analyze,'both')
        colors = [0.4667 0.6745 0.1882];
    end
    
    if strcmpi(condition2analyze,'both')
        % plot cortical surface
        try
            [hFig, iDS, iFig] = view_surface(SurfaceFile);
        catch
            [hFig, iDS, iFig] = view_surface(SurfaceFile);
        end
        % set background color
        hFig.Color = [1 1 1];
        hold on;
        
        % plot statistically-significant electrodes
        sh  = scatter3(datamat(:,1),datamat(:,2),datamat(:,3),60,'filled');
        sh.CData = colors;
        sh.MarkerFaceAlpha = .5;
        sh.MarkerEdgeColor  = [0 0 0];
        sh.MarkerEdgeAlpha = .5;
        sh.SizeData = 50;
    else
        % get clusters
        Clusters = dbscan(datamat,mindist,minpoints);
        tmptable = sortrows(tabulate(Clusters),-3)
        deleteThis   = 0;
        % delete all electrodes that do not belong to the identified cluster or
        % to cluster below the minimum number of electrodes set for the cluster
        % analysis
        if any(tmptable(:,2) < minpoints)
            deleteThis = tmptable(tmptable(:,2) < minpoints,1);
        end
        
        % delete clusters where n_electrodes < minpoints
        datamat(ismember(Clusters,[-1; deleteThis]),:,:)  = [];
        subIDelec(ismember(Clusters,[-1; deleteThis]),:,:)  = [];
        ValRange(ismember(Clusters,[-1; deleteThis]),:,:)  = [];
        Clusters(ismember(Clusters,[-1; deleteThis])) = [];
        save([data_dir,filesep,'dbscan_results_' condition2analyze '_' band2analyze '.mat'],'datamat','ValRange','subIDelec')

        % and combine smaller surviving clusters (where n_electrodes > minpoints) into one
        if size(tmptable,1)-1 > 3
            Clusters(ismember(Clusters,tmptable(4:end,1))) = 3;
        end
        
        % plot cortical surface
        try
            [hFig, iDS, iFig] = view_surface(SurfaceFile);
        catch
            [hFig, iDS, iFig] = view_surface(SurfaceFile);
        end
        % set background color
        hFig.Color = [1 1 1];
        hold on;
        
        % plot statistically-significant electrodes
        clustersID = unique(Clusters);
        sh  = scatter3(datamat(:,1),datamat(:,2),datamat(:,3),60,'filled');
        sh.CData = colors;
        sh.MarkerFaceAlpha = .5;
        sh.MarkerEdgeColor  = [0 0 0];
        sh.MarkerEdgeAlpha = .5;
        sh.SizeData = 50;
        % print parameters used as title
        % title(['mindist = ' num2str(mindist) ', minelecs = ' num2str(minpoints)],'FontWeight','normal','FontSize',22);
        
        h = zeros(1, 1);
        h(1) = plot(NaN,NaN,'color',colors(1,:));
    end
end