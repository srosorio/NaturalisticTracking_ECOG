clear, clc

% music
load('Electrode_AnatLabels_Music_SFB.mat')
Music_SFB = sortrows(tabulate(cellstr(AnatLabels)),1,'ascend');
Music_SFB(2,:) = [];
load('Electrode_AnatLabels_Music_HFB.mat')
Music_HFB = sortrows(tabulate(cellstr(AnatLabels)),1,'ascend');
Music_HFB(4,:) = [];

Values = [[cell2mat(Music_SFB(1,2));0;cell2mat(Music_SFB(2:3,2));0;cell2mat(Music_SFB(4:end,2))], cell2mat(Music_HFB(:,2))];
figure(1),clf,
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
set(gca,'FontSize',12, 'ylim', [0 50])

% speech
load('Electrode_AnatLabels_Speech_SFB.mat')
Speech_SFB = sortrows(tabulate(cellstr(AnatLabels)),1,'ascend');
load('Electrode_AnatLabels_Speech_HFB.mat')
Speech_HFB = sortrows(tabulate(cellstr(AnatLabels)),1,'ascend');

Values = [cell2mat(Speech_SFB(:,2)), [cell2mat(Speech_HFB(1:3,2));0;cell2mat(Speech_HFB(4:end,2))]];
figure(2),clf,
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
set(gca,'FontSize',12, 'ylim', [0 250])