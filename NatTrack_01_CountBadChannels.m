clear, clc,

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
            
for sub_i=1:length(sub2plot)
    % get data from brainstorm database
    cd(['F:\Matlab\brainstorm_db\iEEG\data\' sub2plot{sub_i}]);
    FileTask = dir([sub2plot{sub_i} '*film*']);

    % make sure to keep only good channels
    load([cd,filesep,FileTask.name,filesep,'data_block001.mat'], 'ChannelFlag');
    load([cd,filesep,FileTask.name,filesep,'channel.mat'], 'Channel');
    
    ChannelLabelTask  = squeeze(struct2cell(Channel))';
    Channels2keepTask = strcmp(ChannelLabelTask(:,3),'ECOG');       % keep only ECOG channels
    ChannelFlagTask   = ChannelFlag(Channels2keepTask);
    n_bads(sub_i) = sum(ChannelFlagTask == -1);
    
    clear ChannelFlag Channel
    
end