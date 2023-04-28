% cd('E:\Matlab\brainstorm3');
% brainstorm
% cd('E:\Matlab\IEEG\Scripts');
%%
clear, clc,
segestimation     = 1;
condition2analyze = 'music';
band2analyze      = 'HFB';

if strcmpi(condition2analyze,'music')
    colors = [0.7176 0.2745 1.0000; ...
              0.4902 0.2706 1.0000; ...
              0.2706 0.3686 1.0000];
elseif strcmpi(condition2analyze,'speech')
    colors = [1.0000 0.4118 0.1608; ...
              1.0000 0.6784 0.1608; ...
              0.9412 0.8392 0.3373]; 
elseif strcmpi(condition2analyze,'both')
    colors = [0.4667 0.6745 0.1882; ...
              0.1882 0.6706 0.5255; ...
              0.1882 0.5569 0.6706];    
end
    
% load data to plot
if segestimation == 1
    load(['E:\Matlab\IEEG\Data\CROSdata_',band2analyze,'_windowed.mat']);
else
    load(['E:\Matlab\IEEG\Data\CROSdata_',band2analyze,'.mat']);
end    

% locate data in separate variables
rhos4speech = dataMat(:,:,1);
subeffect   = sum(any(rhos4speech));
disp(['SPEECH: Statistically significant data in ' num2str(subeffect) ' out of ' num2str(length(sub2plot)) ' subjects'])

rhos4music  = dataMat(:,:,2);
subeffect   = sum(any(rhos4music));
disp(['MUSIC: Statistically significant data in ' num2str(subeffect) ' out of ' num2str(length(sub2plot)) ' subjects'])

% Plot MNI surface using brainstorm
SurfaceFile   = 'E:\MATLAB\brainstorm_db\iEEG\anat\@default_subject\tess_cortex_pial_low.mat'; %['C:\Users\andre\OneDrive\Documentos\MATLAB\brainstorm_db\iEEG\anat\' sub2plot '\tess_cortex_central_low.mat'];
%%
if strcmpi(condition2analyze,'speech')
    counter = 1;
    
    for sub_i=1:length(sub2plot)
        
        clear ThisSubStruct
        ThisSubStruct = cell2struct(AllChannelLabels{sub_i},names4fields,2);
        
        for idx=1:length(dataMat)
            if ~isnan(rhos4speech(idx,sub_i)) %&& isnan(rhos4music(idx,sub_i)) %%
                %             test(sub_i,idx,:)    = [ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3)];
                testMat(counter,:)   = [ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3)];
                ValRange(counter,:)  = rhos4speech(idx,sub_i);
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
            if  ~isnan(rhos4music(idx,sub_i)) %% isnan(rhos4speech(idx,sub_i)) &&
                %             test(sub_i,idx,:)    = [ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3)];
                testMat(counter,:)   = [ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3)];
                ValRange(counter,:)  = rhos4speech(idx,sub_i);
                subIDelec(counter,:) = [sub_i,idx];
                counter = counter + 1;
            end
        end
    end  
elseif strcmpi(condition2analyze,'both')
    counter = 1;    
    
    for sub_i=1:length(sub2plot)
        
        clear ThisSubStruct
        ThisSubStruct = cell2struct(AllChannelLabels{sub_i},names4fields,2);
        
        for idx=1:length(dataMat)
            if  ~isnan(rhos4music(idx,sub_i)) && ~isnan(rhos4speech(idx,sub_i)) 
                %             test(sub_i,idx,:)    = [ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3)];
                testMat(counter,:)   = [ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3)];
                ValRange(counter,:)  = rhos4speech(idx,sub_i);
                subIDelec(counter,:) = [sub_i,idx];
                counter = counter + 1;
            end
        end
    end  
end

% close
[hFig, iDS, iFig] = view_surface(SurfaceFile);
hFig.Color = [1 1 1];
hold on;

sh  = scatter3(testMat(:,1),testMat(:,2),testMat(:,3),60,'filled');
sh.MarkerFaceAlpha = .8;
sh.SizeData = 100;

%%
clear testClusters prctelecs numclusters numoutliers
dsts   = 0:0.002:0.04;
points = 12:2:20;

for idx=1:length(dsts)
    clear clusterdata
    for jdx=1:length(points)
        tmp = dbscan(testMat,dsts(idx),points(jdx));
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


figure(3), clf
plot(dsts,(prctelecs ./ (dsts*100) - (prctelecs ./ points')),'LineWidth',1.5)
ylabel('electrodes (%)');
xlabel('\epsilon');
legend(num2str(points'));
ylim([0 100]); set(gca,'FontSize',22);
if strcmpi(band2analyze,'SFB')
    band4title = '1-8 hz';
elseif strcmpi(band2analyze,'HFB')
    band4title = '70-120 hz';
end
title([condition2analyze ' - ' band4title],'FontWeight','normal')
box off
%%
mindist   = 0.02;
minpoints = 12;

testClusters = dbscan(testMat,mindist,minpoints);
tmptable = tabulate(testClusters)
deleteThis   = 0;

if any(tmptable(:,2) < minpoints)
    deleteThis = tmptable(find(tmptable(:,2) < minpoints),1);
end

testMat(testClusters == -1 | testClusters == deleteThis,:,:)   = [];
testClusters(testClusters == -1 | testClusters == deleteThis) = []; 

% close
[hFig, iDS, iFig] = view_surface(SurfaceFile);
hFig.Color = [1 1 1];
hold on;


clustersID = unique(testClusters);

for cluster_i = 1:length(clustersID)
    ThisCluster = testMat(testClusters == clustersID(cluster_i),:);
    sh  = scatter3(ThisCluster(:,1),ThisCluster(:,2),ThisCluster(:,3),60,'filled');
    sh.CData = colors(cluster_i,:); 
    sh.MarkerFaceAlpha = .8;
    sh.SizeData = 100;
end

%%
% if strcmpi(condition2analyze,'both')
%     title(['music and speech - ' band4title],'FontWeight','normal','FontSize',22);   
% else
%     title([condition2analyze ' - ' band4title],'FontWeight','normal','FontSize',22);
% end
% 
% h = zeros(1, 1);
% h(1) = plot(NaN,NaN,'color',colors(1,:));
% h(2) = plot(NaN,NaN,'color',colors(2,:));
% h(3) = plot(NaN,NaN,'color',colors(3,:));
% 
% if length(clustersID) == 1
%     [lh,icons] = legend(h, 'cluster 1','box','off','FontSize',24);
%     icons(2).LineWidth = 5;
% elseif length(clustersID) == 2
%     [lh,icons] = legend(h, 'cluster 1','cluster 2','box','off','FontSize',24);
%     icons(3).LineWidth = 5;
%     icons(5).LineWidth = 5;
% elseif length(clustersID) == 3
%     [lh,icons] = legend(h, 'cluster 1','cluster 2','cluster 3','box','off','FontSize',24);
%     icons(4).LineWidth = 5;
%     icons(6).LineWidth = 5;
%     icons(8).LineWidth = 5;    
% end
% %%
% 
% [hFig, iDS, iFig] = view_surface(SurfaceFile);
% hFig.Color = [1 1 1];
% hold on;
% 
% for sub_i=1:length(sub2plot)
%    
%     clear ThisSubStruct
%     ThisSubStruct = cell2struct(AllChannelLabels{sub_i},names4fields,2);
%     
%     for idx=1:length(dataMat) %sum(~isnan(dataMat(:,sub_i,1)))
%         if ~isnan(rhos4speech(idx,sub_i)) && isnan(normDTWmusic(idx,sub_i)) %%
%             sh  = scatter3(ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3),60,'filled');
% %             txt = ThisSubStruct(idx).Name;
% %             text(ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3),txt);
%             sh.CData = [0.6353    0.0784    0.1843];  % red speech
%             sh.MarkerFaceAlpha = .5;
%             sh.SizeData = rhos4speech(idx,sub_i);
% %             sh.SizeData = 100;
%         elseif isnan(rhos4speech(idx,sub_i)) && ~isnan(normDTWmusic(idx,sub_i))
%             sh = scatter3(ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3),60,'filled');
% %             txt = ThisSubStruct(idx).Name;
% %             text(ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3),txt);
%             sh.CData = [0    0.4471    0.7412];    % blue music
%             sh.MarkerFaceAlpha = .5;
%             sh.SizeData = normDTWmusic(idx,sub_i);
% %             sh.SizeData = 100;
%         elseif ~isnan(rhos4speech(idx,sub_i)) && ~isnan(normDTWmusic(idx,sub_i)) 
%             meanVal = mean([rhos4speech(idx,sub_i), normDTWmusic(idx,sub_i)]);
%             sh = scatter3(ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3),60,'filled');
% %             txt = ThisSubStruct(idx).Name;
% %             text(ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3),txt);
%             sh.CData = [0.9294    0.6941    0.1255];      % green both
%             sh.MarkerFaceAlpha = .5;
%             sh.SizeData = meanVal;      
% %             sh.SizeData = 100;
%         end
%     end
% end
% %%
% % display % of electrodes that survive statistics
% plottedElecs = sum([sum(~isnan(rhos4speech),'all') sum(~isnan(normDTWmusic),'all')]);
% disp([num2str((plottedElecs*100)/TotalElecs) '% electrodes out of ' num2str(TotalElecs) ' show a significant effect']);
% 
% h = zeros(3, 1);
% h(1) = plot(NaN,NaN,'color',[0.6353    0.0784    0.1843]);
% h(2) = plot(NaN,NaN,'color', [0    0.4471    0.7412]);
% h(3) = plot(NaN,NaN,'color',[0.9294    0.6941    0.1255]);
% [lh,icons] = legend(h, 'speech','music','speech + music','box','off','FontSize',24);
% icons(4).LineWidth = 5; 
% icons(6).LineWidth = 5;
% icons(8).LineWidth = 5;
% if strcmpi(band2analyze,'Theta')
%     th = title('Delta/Theta (1-8 Hz)','FontWeight','Normal','FontSize',28);
% else
%     th = title('HFB (70-120 Hz)','FontWeight','Normal','FontSize',28);
% end    
% th.Position = [-0.15   0.0002    -0.015];
%%

% SurfaceFile     = 'E:\Matlab\brainstorm3\defaults\anatomy\ICBM152\tess_cortex_pial_low.mat'; %['C:\Users\andre\OneDrive\Documentos\MATLAB\brainstorm_db\iEEG\anat\' sub2plot '\tess_cortex_central_low.mat'];
% SurfaceFile     = ['C:\Users\andre\OneDrive\Documentos\MATLAB\brainstorm_db\iEEG\anat\' sub2plot '\tess_cortex_central_low.mat']
% ChannelsFile    = 'C:\Users\andre\OneDrive\Documentos\MATLAB\brainstorm_db\iEEG\data\sub-03/@rawsub-03_ses-iemu_ieeg_sub-03_ses-iemu_task-film_acq-clinical_run-1_ieeg/channel.mat';
% load(ChannelsFile);
% close,
% [hFig, iDS, iFig] = view_surface(SurfaceFile);
% hold on;
% 
% for sub_i=1:length(sub2plot)
%     
%     cd(['E:\MATLAB\brainstorm_db\iEEG\data\' sub2plot{sub_i}]);
%     ThisFile = dir([sub2plot{sub_i} '*film*']);
%     cd(ThisFile.name);
%     load('data_block001.mat', 'ChannelFlag');
%     load('channel.mat');
%     
%     cd(['E:\Matlab\IEEG\' sub2plot{sub_i}]);
%         
%     Channel = Channel(ChannelFlag == 1);
%     
%     for idx=1:sum(~isnan(normDTWdata(:,sub_i,1)))
%         if ~isnan(normDTWmusic(idx,sub_i))
%             sh = scatter3(Channel(idx).Loc(1),Channel(idx).Loc(2),Channel(idx).Loc(3),60,'filled');
%             sh.CData = [normDTWmusic(idx,sub_i) 0    0];
%         end
%     end
% 
% end
% %%
% % [hFig, iDS, iFig] = view_channels(ChannelsFile,'ECOG',1,0);
% 
% % hFig.CurrentObject.MarkerSize = 8;
% % hFig.CurrentObject.MarkerEdgeColor = [0 1 0];
% % hFig.CurrentObject.MarkerFaceColor = [0 1 0];
% % set marker colors
% % hFig.CurrentObject.FaceVertexCData = TestColor;
% 
% 
% histogram(normDTWmusic(normDTWmusic ~= 1),10); hold on
% histogram(normDTWspeech(normDTWspeech ~= 1),10)
% legend({'music','speech'},'location','northeast')