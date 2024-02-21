% set paths
clear,clc
iEEG_dir = 'F:\Matlab\IEEG';
data_dir = [iEEG_dir,filesep,'Data'];


tbl = readtable([data_dir,filesep,'AnovaNatTrackTable_wn_ANAT_Music.csv']);
% tbl.lag = tbl.lag / 200;
means = [mean(tbl.r(contains(tbl.frequency,'SFB'))),mean(tbl.r(contains(tbl.frequency,'HFB')))];
err   = [std(tbl.r(contains(tbl.frequency,'SFB')))/sqrt(length(tbl.r(contains(tbl.frequency,'SFB')))),std(tbl.r(contains(tbl.frequency,'HFB')))/sqrt(length(tbl.r(contains(tbl.frequency,'SFB'))))];

clf,subplot(1,3,1)
bh = bar(means, 'FaceColor','Flat','BarWidth', .5); 
hold on;
% bh.FaceColor = [0.8157 0.6000 0.9490];
bh.EdgeColor = [1 1 1];
bh.CData(1,:) = [0.7176 0.2745 1.0000];
bh.CData(2,:) = [0.8157 0.6000 0.9490];
bh.FaceAlpha = .8;
ylim([0 .16]); xlim([0 3]); yticks(0:0.04:0.16)
xtickangle(45)
errorbar(1,means(1),err(1),'k')
errorbar(2,means(2),err(2),'k')
xticklabels({'SFB','HFB'})
title('Music', 'FontWeight', 'Normal');
ylabel('Mean cross-correlation (r)')
set(gca,'FontSize',12)
box off
axis square

clear, clc
iEEG_dir = 'F:\Matlab\IEEG';
data_dir = [iEEG_dir,filesep,'Data'];

tbl = readtable([data_dir,filesep,'AnovaNatTrackTable_wn_ANAT_Speech.csv']);
% tbl.lag = tbl.lag / 200;
means = [[mean(tbl.r(contains(tbl.frequency,'SFB') & contains(tbl.region,'ifg'))) ; ...
         mean(tbl.r(contains(tbl.frequency,'HFB') & contains(tbl.region,'ifg')))], ...
         [mean(tbl.r(contains(tbl.frequency,'SFB') & contains(tbl.region,'mtg'))); ...
         mean(tbl.r(contains(tbl.frequency,'HFB') & contains(tbl.region,'mtg')))], ...
         [mean(tbl.r(contains(tbl.frequency,'SFB') & contains(tbl.region,'somatomotor'))); ...
         mean(tbl.r(contains(tbl.frequency,'HFB') & contains(tbl.region,'somatomotor')))], ...
         [mean(tbl.r(contains(tbl.frequency,'SFB') & contains(tbl.region,'somatosensory'))); ...
         mean(tbl.r(contains(tbl.frequency,'HFB') & contains(tbl.region,'somatosensory')))], ...
         [mean(tbl.r(contains(tbl.frequency,'SFB') & contains(tbl.region,'stg'))); ...
         mean(tbl.r(contains(tbl.frequency,'HFB') & contains(tbl.region,'stg')))], ...
         [mean(tbl.r(contains(tbl.frequency,'SFB') & contains(tbl.region,'supramarginal'))); ...
         mean(tbl.r(contains(tbl.frequency,'HFB') & contains(tbl.region,'supramarginal')))]];

     
err = [[std(tbl.r(contains(tbl.frequency,'SFB') & contains(tbl.region,'ifg'))) / sqrt(length(tbl.r(contains(tbl.frequency,'SFB') & contains(tbl.region,'ifg')))); ...
         std(tbl.r(contains(tbl.frequency,'HFB') & contains(tbl.region,'ifg'))) / sqrt(length(tbl.r(contains(tbl.frequency,'HFB') & contains(tbl.region,'ifg'))))], ...
         [std(tbl.r(contains(tbl.frequency,'SFB') & contains(tbl.region,'mtg'))) / sqrt(length(tbl.r(contains(tbl.frequency,'SFB') & contains(tbl.region,'mtg')))); ...
         std(tbl.r(contains(tbl.frequency,'HFB') & contains(tbl.region,'mtg'))) / sqrt(length(tbl.r(contains(tbl.frequency,'HFB') & contains(tbl.region,'mtg'))))], ...
         [std(tbl.r(contains(tbl.frequency,'SFB') & contains(tbl.region,'somatomotor'))) / sqrt(length(tbl.r(contains(tbl.frequency,'SFB') & contains(tbl.region,'somatomotor')))); ...
         std(tbl.r(contains(tbl.frequency,'HFB') & contains(tbl.region,'somatomotor'))) / sqrt(length(tbl.r(contains(tbl.frequency,'HFB') & contains(tbl.region,'somatomotor'))))], ...
         [std(tbl.r(contains(tbl.frequency,'SFB') & contains(tbl.region,'somatosensory'))) / sqrt(length(tbl.r(contains(tbl.frequency,'SFB') & contains(tbl.region,'somatosensory')))); ...
         std(tbl.r(contains(tbl.frequency,'HFB') & contains(tbl.region,'somatosensory'))) / sqrt(length(tbl.r(contains(tbl.frequency,'HFB') & contains(tbl.region,'somatosensory'))))], ...
         [std(tbl.r(contains(tbl.frequency,'SFB') & contains(tbl.region,'stg'))) / sqrt(length(tbl.r(contains(tbl.frequency,'SFB') & contains(tbl.region,'stg')))); ...
         std(tbl.r(contains(tbl.frequency,'HFB') & contains(tbl.region,'stg')))] / sqrt(length(tbl.r(contains(tbl.frequency,'HFB') & contains(tbl.region,'stg')))), ...
         [std(tbl.r(contains(tbl.frequency,'SFB') & contains(tbl.region,'supramarginal'))) / sqrt(length(tbl.r(contains(tbl.frequency,'SFB') & contains(tbl.region,'supramarginal')))); ...
         std(tbl.r(contains(tbl.frequency,'HFB') & contains(tbl.region,'supramarginal')))]] / sqrt(length(tbl.r(contains(tbl.frequency,'HFB') & contains(tbl.region,'supramarginal'))));

labels = {'IFG','MTG', 'Somatomotor', 'Somatosensory', 'STG', 'Supramarginal'};

subplot(1,3,2)
bh = bar(means', 'FaceColor', 'Flat') ;    
xticklabels(labels)
ylabel('Mean cross-correlation (r)');
xtickangle(45)
bh(1).EdgeColor = [1 1 1];
bh(1).CData = repmat([1.0000 0.4118 0.1608],6,1);
bh(2).EdgeColor = [1 1 1];
bh(2).CData = repmat([0.9804 0.5804 0.4118],6,1);
bh(1).FaceAlpha = .8;
bh(2).FaceAlpha = .8;
legend({'SFB','HFB'},'box', 'off', 'Orientation', 'horizontal', 'Location', 'north','AutoUpdate','off')
title('Speech', 'FontWeight', 'Normal');
ylim([0 .16]); xlim([0 7]); 
yticks(0:0.04:0.16)
box off
hold on
set(gca,'FontSize',12)

for k = 1:numel(bh)                                                      % Earlier MATLAB Versions
    ctr(k,:) = bsxfun(@plus, bh(k).XData, [bh(k).XOffset]');
    ydt(k,:) = bh(k).YData;
end

for idx=1:length(ctr)
    errorbar(ctr(1,idx), ydt(1,idx), err(1,idx), 'k', 'MarkerSize',0.1)
    errorbar(ctr(2,idx), ydt(2,idx), err(2,idx), 'k', 'MarkerSize',0.1)
end

bh(1).BarWidth = 1;
bh(2).BarWidth = 1;
axis square

%
means = [[mean(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.region,'ifg'))) ; ...
         mean(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.region,'ifg')))], ...
         [mean(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.region,'mtg'))); ...
         mean(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.region,'mtg')))], ...
         [mean(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.region,'somatomotor'))); ...
         mean(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.region,'somatomotor')))], ...
         [mean(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.region,'somatosensory'))); ...
         mean(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.region,'somatosensory')))], ...
         [mean(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.region,'stg'))); ...
         mean(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.region,'stg')))], ...
         [mean(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.region,'supramarginal'))); ...
         mean(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.region,'supramarginal')))]];

err = [[std(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.region,'ifg'))) / sqrt(length(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.region,'ifg')))); ...
         std(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.region,'ifg'))) / sqrt(length(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.region,'ifg'))))], ...
         [std(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.region,'mtg'))) / sqrt(length(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.region,'mtg')))); ...
         std(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.region,'mtg'))) / sqrt(length(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.region,'mtg'))))], ...
         [std(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.region,'somatomotor'))) / sqrt(length(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.region,'somatomotor')))); ...
         std(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.region,'somatomotor'))) / sqrt(length(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.region,'somatomotor'))))], ...
         [std(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.region,'somatosensory'))) / sqrt(length(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.region,'somatosensory')))); ...
         std(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.region,'somatosensory'))) / sqrt(length(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.region,'somatosensory'))))], ...
         [std(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.region,'stg'))) / sqrt(length(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.region,'stg')))); ...
         std(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.region,'stg')))] / sqrt(length(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.region,'stg')))), ...
         [std(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.region,'supramarginal'))) / sqrt(length(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.region,'supramarginal')))); ...
         std(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.region,'supramarginal')))]] / sqrt(length(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.region,'supramarginal'))));

     labels = {'IFG','MTG', 'Somatomotor', 'Somatosensory', 'STG', 'Supramarginal'};

subplot(1,3,3)
bh = bar(means', 'FaceColor', 'Flat') ;    
xticklabels(labels)
ylabel('Mean lag (s)');
xtickangle(45)
bh(1).EdgeColor = [1 1 1];
bh(1).CData = repmat([1.0000 0.4118 0.1608],6,1);
bh(2).EdgeColor = [1 1 1];
bh(2).CData = repmat([0.9804 0.5804 0.4118],6,1);
bh(1).FaceAlpha = .8;
bh(2).FaceAlpha = .8;
legend({'SFB','HFB'},'box', 'off', 'Orientation', 'horizontal', 'Location', 'north','AutoUpdate','off')
title('Speech', 'FontWeight', 'Normal');
ylim([-.2 1]); xlim([0 7]); 
% yticks(0:0.04:0.16)
box off
hold on
plot(get(gca,'xlim'),[0 0],'k');
set(gca,'FontSize',12)

for k = 1:numel(bh)                                                      % Earlier MATLAB Versions
    ctr(k,:) = bsxfun(@plus, bh(k).XData, [bh(k).XOffset]');
    ydt(k,:) = bh(k).YData;
end

for idx=1:length(ctr)
    errorbar(ctr(1,idx), ydt(1,idx), err(1,idx), 'k', 'MarkerSize',0.1)
    errorbar(ctr(2,idx), ydt(2,idx), err(2,idx), 'k', 'MarkerSize',0.1)
end

bh(1).BarWidth = 1;
bh(2).BarWidth = 1;
axis square
