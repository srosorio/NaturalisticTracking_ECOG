clear, clc
iEEG_dir = 'F:\Matlab\IEEG';
data_dir = [iEEG_dir,filesep,'Data'];

% music
load([data_dir,filesep,'Electrode_AnatLabels_Music_SFB.mat'])
Music_SFB = sortrows(tabulate(cellstr(AnatLabels)),1,'ascend');
Music_SFB(2,:) = [];
load([data_dir,filesep,'Electrode_AnatLabels_Music_HFB.mat'])
Music_HFB = sortrows(tabulate(cellstr(AnatLabels)),1,'ascend');
Music_HFB(4,:) = [];

Values = [[cell2mat(Music_SFB(1,2));0;cell2mat(Music_SFB(2:3,2));0;cell2mat(Music_SFB(4:end,2))], cell2mat(Music_HFB(:,2))];
subplot(1,3,2)
bh = bar(Values,'BarWidth',1,'FaceColor','flat');
bh(1).CData = repmat([0.7176 0.2745 1.0000],9,1);
bh(2).CData = repmat([0.8157 0.6000 0.9490],9,1);
bh(1).FaceAlpha = .8; bh(2).FaceAlpha = .8; 
bh(1).EdgeColor = [1 1 1]; bh(2).EdgeColor = [1 1 1];
xticklabels({'IFG','IPC','MTG','PFC','Premotor','Somatomotor','Somatosensory','STG','Supramarginal'})
xtickangle(45)
legend({'SFB','HFB'}, 'box', 'off', 'Orientation', 'horizontal', 'Location', 'north')
ylabel('Electrode count')
title('Music','FontWeight','normal')
box off
set(gca,'FontSize',12, 'ylim', [0 50]);
axis square
tmptbl = table({'IFG','IPC','MTG','PFC','Premotor','Somatomotor','Somatosensory','STG','Supramarginal'}', ...
    Values(:,1), Values(:,2), 'VariableNames',{'Region','SFB','HFB'});
writetable(tmptbl,[data_dir,filesep,'ElecCountPerAnatRegion_Music.xlsx']);

% speech
load([data_dir,filesep,'Electrode_AnatLabels_Speech_SFB.mat'])
Speech_SFB = sortrows(tabulate(cellstr(AnatLabels)),1,'ascend');
load([data_dir,filesep,'Electrode_AnatLabels_Speech_HFB.mat'])
Speech_HFB = sortrows(tabulate(cellstr(AnatLabels)),1,'ascend');

Values = [cell2mat(Speech_SFB(:,2)), [cell2mat(Speech_HFB(1:3,2));0;cell2mat(Speech_HFB(4:end,2))]];
subplot(1,3,3)
bh = bar(Values,'BarWidth',1,'FaceColor','flat');
bh(1).CData = repmat([1.0000 0.4118 0.1608],8,1);
bh(2).CData = repmat([0.9804 0.5804 0.4118],8,1);
bh(1).FaceAlpha = .8; bh(2).FaceAlpha = .8; 
bh(1).EdgeColor = [1 1 1]; bh(2).EdgeColor = [1 1 1];
xticklabels({'IFG','MTG','PFC','Premotor','Somatomotor','Somatosensory','STG','Supramarginal'})
xtickangle(45)
legend({'SFB','HFB'}, 'box', 'off', 'Orientation', 'horizontal', 'Location', 'north')
ylabel('Electrode count')
title('Speech','FontWeight','normal')
box off
set(gca,'FontSize',12, 'ylim', [0 250], 'xlim',[0 9])
axis square

tmptbl = table({'IFG','MTG','PFC','Premotor','Somatomotor','Somatosensory','STG','Supramarginal'}' ...
    ,Values(:,1), Values(:,2), 'VariableNames',{'Region','SFB','HFB'});
writetable(tmptbl,[data_dir,filesep,'ElecCountPerAnatRegion_Speech.xlsx']);
%%
% both
clf
subplot(1,3,1);
Music_SFB = Music_SFB(contains(Music_SFB(:,1),{'stg','mtg','supramarginal'}),:);
Music_HFB = Music_HFB(contains(Music_HFB(:,1),{'stg','mtg','supramarginal'}),:);
Speech_SFB = Speech_SFB(contains(Speech_SFB(:,1),{'stg','mtg','supramarginal'}),:);
Speech_HFB = Speech_HFB(contains(Speech_HFB(:,1),{'stg','mtg','supramarginal'}),:);

Values = [cell2mat(Music_SFB(:,2)), cell2mat(Music_HFB(:,2)); cell2mat(Speech_SFB(:,2)), cell2mat(Speech_HFB(:,2)) ];
bh = bar(Values,'BarWidth',1,'FaceColor','flat');
xticklabels({'STG','MTG','Supramarginal','STG','MTG','Supramarginal'})
xtickangle(45)
% legend({'SFB','HFB'}, 'box', 'off', 'Orientation', 'horizontal', 'Location', 'north')
set(gca,'FontSize',12, 'ylim', [0 250], 'xlim',[0 7])
bh(1).CData = [repmat([0.7176 0.2745 1.0000],3,1); repmat([1.0000 0.4118 0.1608],3,1)];
bh(2).CData = [repmat([0.8157 0.6000 0.9490],3,1); repmat([ 0.9804 0.5804 0.4118],3,1)];
bh(1).FaceAlpha = .8; bh(2).FaceAlpha = .8; 
bh(1).EdgeColor = [1 1 1]; bh(2).EdgeColor = [1 1 1];
ylabel('Electrode count');
box off
axis square

tmptbl = table({'STG','MTG','Supramarginal','STG','MTG','Supramarginal'}', ...
    Values(:,1), Values(:,2), 'VariableNames',{'Region','SFB','HFB'});
writetable(tmptbl,[data_dir,filesep,'ElecCountPerAnatRegion_Both.xlsx']);

tbl = readtable([data_dir,filesep,'AnovaNatTrackTable_wn_ANAT.csv']);
means = [mean(tbl.r(contains(tbl.condition,'Music'))),mean(tbl.r(contains(tbl.condition,'Speech')))];
err   = [std(tbl.r(contains(tbl.condition,'Music')))/sqrt(length(tbl.r(contains(tbl.condition,'Music')))),std(tbl.r(contains(tbl.condition,'Speech')))/sqrt(length(tbl.r(contains(tbl.condition,'Speech'))))];

subplot(1,3,2)
bh = bar(means, 'FaceColor','Flat','BarWidth', .5); 
hold on;
bh.EdgeColor = [1 1 1];
bh.CData(1,:) = [0.7176 0.2745 1.0000];
bh.CData(2,:) = [1.0000 0.4118 0.1608];
bh.FaceAlpha = .8;
ylim([0 .16]); xlim([0 3]); yticks(0:0.04:0.16)
xtickangle(45)
errorbar(1,means(1),err(1),'k')
errorbar(2,means(2),err(2),'k')
xticklabels({'Music','Speech'})
ylabel('Mean cross-correlation (r)')
set(gca,'FontSize',12)
box off
axis square
%

means = [[mean(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.condition,'Music'))) ; ...
         mean(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.condition,'Music')))], ...
         [mean(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.condition,'Speech'))); ...
         mean(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.condition,'Speech')))]];
     
err = [[std(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.condition,'Music'))) / sqrt(length(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.condition,'Music')))); ...
         std(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.condition,'Music'))) / sqrt(length(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.condition,'Music'))))], ...
         [std(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.condition,'Speech'))) / sqrt(length(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.condition,'Speech')))); ...
         std(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.condition,'Speech'))) / sqrt(length(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.condition,'Speech'))))]];

subplot(1,3,3)
bh = bar(means, 'FaceColor','Flat','BarWidth', .5); 
hold on;
bh(1).EdgeColor = [1 1 1];
bh(2).EdgeColor = [1 1 1];
bh(1).CData = [0.7176 0.2745 1.0000; 1.0000 0.4118 0.1608];
bh(2).CData = [0.8157 0.6000 0.9490; 0.9804 0.5804 0.4118];
bh(1).FaceAlpha = .8;
bh(2).FaceAlpha = .8;
ylim([-.2 1]); 
xtickangle(45)
plot(get(gca,'xlim'),[0 0],'k');

for k = 1:numel(bh)                                                      % Earlier MATLAB Versions
    ctr(k,:) = bsxfun(@plus, bh(k).XData, [bh(k).XOffset]');
    ydt(k,:) = bh(k).YData;
end

for idx=1:length(ctr)
    errorbar(ctr(1,idx), ydt(1,idx), err(1,idx), 'k', 'MarkerSize',0.1)
    errorbar(ctr(2,idx), ydt(2,idx), err(2,idx), 'k', 'MarkerSize',0.1)
end
xticklabels({'Music','Speech'})
ylabel('Mean lag (s)')
set(gca,'FontSize',12)
box off
axis square