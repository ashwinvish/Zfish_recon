function [len] = treeLength(cellID)
%% Calculate tree length of tree associated with cellID

% file path to location of .swc files
filePath  = '/Users/ashwin/Google Drive/Zfish/LowEMtoHighEM//SWC_all/combinedConsensus-resampled/';
fileName = sprintf('%d_reRoot_reSample_5000.swc',cellID);

% load_tree, len_tree funcitons are from the treestoolbox
% (https://www.treestoolbox.org/)
[inTree,~,~] = load_tree(fullfile(filePath,fileName));
len = [0;len_tree(inTree)];
len = sum(len)/1000; % covert to microns
%disp('no tree')

end