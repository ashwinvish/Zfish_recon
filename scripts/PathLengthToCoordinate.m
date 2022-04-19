function [pathLength,k] = PathLengthToCoordinate(pointCloud,cellID)
%% Calculate pathlength from root node to coordinates in pointCloud
% pointCloud is mx3 vector
% cellID is the ID of the cell

% file path to .swc files
filePath  = '/Users/ashwin/Google Drive/Zfish/LowEMtoHighEM//SWC_all/combinedConsensus-resampled/';
fileName = sprintf('%d_reRoot_reSample_5000.swc',cellID);
%fileName = sprintf('%d.swc',cellID);

if exist(fullfile(filePath,fileName))==0
    filePath  = '/Users/ashwin/Google Drive/Zfish/LowEMtoHighEM/SWC_all/consensus-20190923/swc/';
end

[inTree,~,~] = load_tree(fullfile(filePath,fileName));

B = [inTree.X,inTree.Y,inTree.Z]./[5,5,45];
len = [0;len_tree(inTree)];
Idpar = ipar_tree(inTree);
pathLength = zeros(size(pointCloud,1),1);
pvec = Pvec_tree(inTree);

k = dsearchn(B,pointCloud);
% parfor i = 1:size(pointCloud,1)
%     [k(i),~] = dsearchn(B,pointCloud(i,:));
%     %pathLength(i) = sum(len(Idpar(k,:)+1));
%     pathLength(i) = pvec(k(i));
% end

pathLength = zeros(length(k),1);

for i = 1:length(k)
    if k(i)~=1
      pathLength(i) = pvec(k(i));
    else
      pathLength(i) = pvec(2);
    end
end
%pathLength = pvec(k);
pathLength = pathLength./1000; % convert to microns
end
