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
plot_neteffect    = 0;  % plot effect prior to statistics
dbscan_param      = 1;  % plot dbscan optimization process
plot_clusters     = 1;  % only after dbscan optimization
condition2analyze = 'speech';
band2analyze      = 'SFB';
perm_type         = 'wn';   % wn = whitenoise, ts = trialshuffling
% -------------------------------------------------------------------------

% set paths
iEEG_dir = 'F:\Matlab\IEEG';
data_dir = [iEEG_dir,filesep,'Data'];

% colors per condition (1,:) music, (2,:) speech, (3,:) music and speech
colors  = [0.7176 0.2745 1.0000; ...
           0.9412 0.5804 0.3373; ...
           0.4667 0.6745 0.1882];              
    
% load data to plot
if strcmpi(perm_type,'ts')
    load(['F:\Matlab\IEEG\Data\CROSdata_',band2analyze,'_trialshuffle.mat']);
elseif strcmpi(perm_type,'wn')
    load(['F:\Matlab\IEEG\Data\CROSdata_',band2analyze,'_whitenoise.mat']);
end    

% locate data in separate variables and plot
rhos4speech = dataMat(:,:,1);
rhos4music  = dataMat(:,:,2);

% get MNI cortical surface using brainstorm. We will plot data here.
SurfaceFile   = 'F:\MATLAB\brainstorm_db\iEEG\anat\@default_subject\tess_cortex_pial_low.mat'; %['C:\Users\andre\OneDrive\Documentos\MATLAB\brainstorm_db\iEEG\anat\' sub2plot '\tess_cortex_central_low.mat'];

try
    [hFig, iDS, iFig] = view_surface(SurfaceFile);
catch
    [hFig, iDS, iFig] = view_surface(SurfaceFile);
end
hFig.Color = [1 1 1];
hold on;
    
n_subs  = length(sub2plot);
n_elecs = length(dataMat);

for sub_i=1:n_subs
    clear ThisSubStruct
    ThisSubStruct = cell2struct(AllChannelLabels{sub_i},names4fields,2);
    n_elecs = length({ThisSubStruct.Name});
    for elec_i=1:n_elecs
        datamat(elec_i,:)   = [ThisSubStruct(elec_i).Loc(1),ThisSubStruct(elec_i).Loc(2),ThisSubStruct(elec_i).Loc(3)];
    end
    elecspersub(sub_i) = length(datamat);
    % plot electrodes
    sh  = scatter3(datamat(:,1),datamat(:,2),datamat(:,3),60,'filled');
    sh.MarkerFaceAlpha = .5;
    sh.SizeData = 50;
    sh.CData    = [0    0.4470    0.7410];
    sh.MarkerEdgeColor  = [0 0 0];
    sh.MarkerEdgeAlpha = .5;
    clear datamat
end

