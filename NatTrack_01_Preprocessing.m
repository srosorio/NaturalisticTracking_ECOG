clear, clc,

addpath('F:\Matlab\fieldtrip-20210709');
addpath('F:\Matlab\IEEG\Scripts\NSL');
ft_defaults;

% create a directory to save data
iEEG_dir = 'F:\Matlab\IEEG';
data_dir = [iEEG_dir,filesep,'Data'];

if ~isfolder(data_dir)
    mkdir(data_dir)
end

% subjects
sub2plot = {'sub-02','sub-03','sub-05','sub-06','sub-10','sub-12','sub-16','sub-18', ...
            'sub-19','sub-20','sub-22','sub-24','sub-25','sub-26','sub-27','sub-34', ...
            'sub-36','sub-36HD','sub-39','sub-40','sub-45','sub-45HD','sub-46','sub-48', ...
            'sub-51','sub-54','sub-55','sub-58','sub-59','sub-60','sub-61','sub-63'};
            
% frequency band to analyze (SLB or HFB)
band2analyze = 'SFB';  
lowfreqs     = [1 8];

% initialize a couple of cell arrays to store data
AllDataStructuresFT = cell(length(sub2plot),2);
AllChannelLabels    = cell(length(sub2plot),1);

for sub_i=1:length(sub2plot)
    % get data from brainstorm database
    cd(['F:\Matlab\brainstorm_db\iEEG\data\' sub2plot{sub_i}]);
    FileTask = dir([sub2plot{sub_i} '*film*']);

    % make sure to keep only good channels
    load([cd,filesep,FileTask.name,filesep,'data_block001.mat'], 'ChannelFlag');
    load([cd,filesep,FileTask.name,filesep,'channel.mat'], 'Channel');
    
    ChannelFlagTask   = ChannelFlag;
    ChannelLabelTask  = squeeze(struct2cell(Channel))';
    Channels2keepTask = strcmp(ChannelLabelTask(:,3),'ECOG');       % keep only ECOG channels
    ChannelLabel2keep = ChannelLabelTask(Channels2keepTask,1);      % get the labels for those channels. We'll need them for rest data
    ChannelFlagTask(~Channels2keepTask) = [];    
    names4fields     = fieldnames(Channel);
    % delete channels that won't be included from the ChannelFlag list
    clear ChannelFlag Channel
    
    cd(['F:\Matlab\IEEG\' sub2plot{sub_i}]);
    
    %% Data preprocessing
    cfg                     = [];
    cfg.dataset             = [sub2plot{sub_i} '.eeg'];
    cfg.headerfile          = [sub2plot{sub_i} '.vmrk'];
    cfg.trialfun            = 'ft_trialfun_general';
    cfg.trialdef.eventtype  = 'Stimulus';
    cfg.trialdef.eventvalue = 'speech';
    cfg.trialdef.prestim    = 0;
    cft.trialdef.poststim   = 30;
    SpeechData              = ft_definetrial(cfg);
    
    cfg.trialdef.eventtype  = 'Stimulus';
    cfg.trialdef.eventvalue = 'music';
    MusicData               = ft_definetrial(cfg);
    
    %% Filter data and resample to match resolution of signal envelopes
    cfg = [];
    cfg.dataset = [sub2plot{sub_i} '.eeg'];
    cfg.trial           = 'all';
    cfg.demean          = 'yes';
    cfg.bpfilter        = 'yes';
    cfg.bpfiltord       = 3;
    if strcmpi(band2analyze,'SFB')
        cfg.bpfreq      = [1 8];
    elseif strcmpi(band2analyze,'HFB')
        cfg.bpfreq      = [70 120];
        cfg.hilbert     = 'abs';
    end
    
    ThisData = ft_preprocessing(cfg);
    
    % delete necessary channels and keep only the ECOG channels identified above
    ThisData.trial{1} = ThisData.trial{1}(Channels2keepTask,:);
    ThisData.label    = ThisData.label(Channels2keepTask);
    % and now keep only good channels
    ThisData.trial{1} = ThisData.trial{1}(ChannelFlagTask == 1,:);
    ThisData.label    = ThisData.label(ChannelFlagTask == 1);    
    
    % define trials according to music and speech segments. Create separate files
    Speech_data = ft_redefinetrial(SpeechData,ThisData);
    Music_data  = ft_redefinetrial(MusicData,ThisData);
    
    % resample data to match sampling rate in cochlear envelope
    cfg = [];
    cfg.resamplefs  = 200;
    Speech_data     = ft_resampledata(cfg,Speech_data);
    Music_data      = ft_resampledata(cfg,Music_data);
    
    % delete the last trial for music data (to match number of trials in speech data)
    Music_data.trial(7) = [];
    Music_data.time(7) = [];
    
    AllDataStructuresFT{sub_i,1} = Speech_data;
    AllDataStructuresFT{sub_i,2} = Music_data;
    
    % and now we keep the channel labels and positions for each subject (only from the Task data file) for the plotting script
    TmpChannelLabelData     = ChannelLabelTask(Channels2keepTask,:);
    AllChannelLabels{sub_i} = TmpChannelLabelData(ChannelFlagTask == 1,:);
end

% save fieldtrip structures
save([data_dir,filesep,'fieldtrip_structures_' band2analyze '.mat'], 'AllDataStructuresFT','AllChannelLabels','names4fields');