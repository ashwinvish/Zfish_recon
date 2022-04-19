function[] = ReRoot(cellID)
% re-root trees allows for visual inspection of an swc tree and re-root if
% needed and resamples the tree to decrease number of nodes.

% cellID is the ID of the cells that is being re-rooted


originalFilePath  = '/Users/ashwin/Google Drive/Zfish/LowEMtoHighEM/SWC_all/consensus-20190415/swc/';
newFilePath = '/Users/ashwin/Google Drive/Zfish/LowEMtoHighEM/SWC_all/combinedConsensus-resampled/';
fileName = sprintf('%d.swc',cellID);

% if exist(fullfile(originalFilePath,fileName))==0
%     %filePath  = '/Users/ashwin/Google Drive/Zfish/LowEMtoHighEM/SWC_all/consensus-20190923/swc/';
%     originalFilePath  = '/Users/ashwin/Google Drive/Zfish/LowEMtoHighEM/SWC_all/consensus-20200406/swc/';
%     newFilePath = '/Users/ashwin/Google Drive/Zfish/LowEMtoHighEM/SWC_all/combinedConsensus-resampled/';
% end

reSampleFactor = 5000; % see trees toolbox for notes

if exist(fullfile(originalFilePath,fileName))
    [tree,~,~] = load_tree(fullfile(originalFilePath,fileName));
    tree = resample_tree(tree,reSampleFactor);
             xplore_tree(tree);
             view(90,0);
    prompt = 'Enter the new root node: '
    %lia = find(tree.X >= max(tree.X)-10000,1);
    %pmt = find(tree.Z == min(tree.Z));
    pmt = 1;
    
    
    if ~isempty(pmt)
        istart = input(prompt);
        %istart = pmt;
        close all;
        [tree, order] = redirect_tree (tree, istart);
        newFileName = sprintf('%d_reRoot_reSample_%d.swc',cellID,reSampleFactor);
        swc_tree(tree, fullfile(newFilePath,newFileName));
    else
        close all;
        newFileName = sprintf('%d_ORoot_reSample_%d.swc',cellID,reSampleFactor);
        swc_tree (tree, fullfile(newFilePath,newFileName));
    end
     
end
