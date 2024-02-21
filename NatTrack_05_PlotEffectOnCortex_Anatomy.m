% initialize BST
cd('F:\Matlab\brainstorm3');
brainstorm
%%

clear, clc
condition2analyze = 'music';   % 'speech' or 'music'
band2analyze      = 'SFB';      % SFB (1-8 Hz) or HFB (70-120 Hz) 

% set paths
iEEG_dir = 'F:\Matlab\IEEG';
data_dir = [iEEG_dir,filesep,'Data'];

% load data
load([data_dir,filesep,'dbscan_results_' condition2analyze '_' band2analyze '.mat']);
load([data_dir,filesep,'Electrode_AnatLabels_' condition2analyze '_' band2analyze '.mat'])

% get MNI cortical surface using brainstorm. We will plot data here.
SurfaceFile   = 'F:\MATLAB\brainstorm_db\iEEG\anat\@default_subject\tess_cortex_pial_low.mat';

%%
[hFig, iDS, iFig] = view_surface(SurfaceFile);
hFig.Color = [1 1 1];
hold on;
cortical_areas = {'somatosensory','somatomotor','ifg','supramarginal','stg','mtg'};
for idx=1:length(cortical_areas)
    plot_this = strcmp(cellstr(AnatLabels),cortical_areas{idx});
    sh = scatter3(datamat(plot_this,1),datamat(plot_this,2),datamat(plot_this,3),60,'filled');
end

