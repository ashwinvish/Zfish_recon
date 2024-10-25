function tree =  SwctoZbrian(cellID)
% SwctoZbrian converts the .swc tree from NG space to Z-brain space.
% cellID is a nx1 vector of the cell IDs that need to be transfromed.

% load path where .swc files are stored
if ismac
    addpath(genpath('/Users/ashwin/Documents/LabWork'));
    fname  = '/Users/ashwin/Google Drive/Zfish/LowEMtoHighEM/SWC_all/combinedConsensus-resampled';
    
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    fname =  '/usr/people/ashwinv/seungmount/research/Ashwin/Z-Brain/LowEMtoHighEM/SWC_all/consensus-20180920/swc/';
end


for i = 1:numel(cellID)   
    fileName = sprintf('%d_reRoot_reSample_5000.swc',cellID(i));
    if exist(fullfile(fname,fileName))
        fID = fopen(fullfile(fname,fileName));
        swc{i} = textscan(fID, '%f %f %f %f %f %f %f','HeaderLines',6,'CollectOutput',true);
        fclose(fID);
        %swc{i} = dlmread(fullfile(fname,fileName), ' ',6,0);
    elseif exist(fullfile(fname,sprintf('%d.swc',cellID(i))))
        fileName = sprintf('%d.swc',cellID(i));
        fID = fopen(fullfile(fname,fileName));
        swc{i} = textscan(fID, '%f %f %f %f %f %f %f','HeaderLines',0,'CollectOutput',true);
        fclose(fID);
    else
        tree{i} = [];
        continue;
    end
    
    coord = swc{i}{1,1}(:, 3:5) ;
    
    % compute offset as data in NG is offset from origin
    offset = [920,752,16400] .* [80, 80, 45]; % offset at mip4 number can be obtained from NG .info file
    
    
    coord(:,1) = coord(:,1) - offset(1);
    coord(:,2) = coord(:,2) - offset(2);
    coord(:,3) = coord(:,3) - offset(3);
    
    % adjust voxel size to LM space
    %voxel --> micron
    
   % coord = double(coord);
    coord = coord./ 1000;
    coord(:,1) = coord(:,1)/0.798;
    coord(:,2) = coord(:,2)/0.798;
    coord(:,3) = coord(:,3)/2;
    
    % rotation NG data is 90 to the righ as compared to the z-brain atlas
    % Create rotation matrix
    theta = pi/2; % to rotate 90 counterclockwise
    R = [cos(theta) -sin(theta) 0 0; sin(theta) cos(theta) 0 0; ...
        0 0 1 0; 0 0 0 1];
    rotateTForm = affine3d(R);
    coord = transformPointsForward(rotateTForm, coord);
    
    % translate the image after it is rotate to compensate for the rotation
    % about the top left corner
    coord(:,2) = coord(:,2) + 436;
    
    % perform transformation
    % load matlab transform object to transform from NG space to Z-brain
    % space
    
    if ismac
        load('/Users/ashwin/Google Drive/Zfish/LowEMtoHighEM/tformRough-LowEMtoHighEM-set2-Elavl3-Mnx-SB-set4.mat');
    else
        load('/usr/people/ashwinv/seungmount/research/Ashwin/Z-Brain/LowEMtoHighEM/tformRough-LowEMtoHighEM-set2.mat');
    end
    
    % transform from voxels to microns
    coord(:,1) = coord(:,1) * 0.798 ;
    coord(:,2) = coord(:,2) * 0.798 ;
    coord(:,3) = coord(:,3) * 2 ;
    
    % transform from NG space to Atlas space using a precomputed
    % transfromation matrix.
    swc_new{i} = swc{i}{1,1};
    swc_new{i}(:,3:5) = transformPointsForward(tformRough, coord);
    
    newFileName = sprintf('%d_reRoot_reSample_5000_transfromed.swc',cellID(i));
    dlmwrite(newFileName,swc_new{i},' ');
    tree{i} = load_tree(newFileName);
    delete(newFileName);  
end

end


