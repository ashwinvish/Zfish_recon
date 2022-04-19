function [pointCloudZB] =  TransformPoints(pointCloudNG,mipLevel)
%TransfromPoints transfroms pointcloud from NG space to Z-brian space
% pointCloudNG is the MX3 point cloud in NG space
% mip level is either 0 or 4
% is 0 if the piont clouds are queried from the edge list
% is 4 if the poit clouds are from the swc skeletons.
% pointCloudZB is the poiuntcloud in Z-brain space.

if mipLevel == 0
    coord = pointCloudNG;
    
    % multiply with mip0 resolution to get to nm space
    coord = coord .* [5,5,45];
    
    % convert the mip0 offset to nm space
    offset = [14720, 12032, 16400] .* [5,5,45]; % offset can be obtained from NG .info file.
    
    % subtract offset to line up with origin
    coord(:,1) = coord(:,1) - offset(1);
    coord(:,2) = coord(:,2) - offset(2);
    coord(:,3) = coord(:,3) - offset(3);
    
    % adjust voxel size to LM space
    % voxel --> micron
    
    coord = coord ./ 1000;
    coord(:,1) = coord(:,1)/0.798;
    coord(:,2) = coord(:,2)/0.798;
    coord(:,3) = coord(:,3)/2;
    
    % rotate the above points to get to z-brian orientation. NG space is
    % rotated 90 to the right as compared to z-brian atlas.
    
    % Create rotation matrix
    theta = pi/2; % to rotate 90 counterclockwise
    R = [cos(theta) -sin(theta) 0 0; sin(theta) cos(theta) 0 0; ...
        0 0 1 0; 0 0 0 1];
    rotateTForm = affine3d(R);
    coord = transformPointsForward(rotateTForm, coord);
    
    % following rotation a translation needs to be applied to bring back to
    % origin.
    
    coord(:,2) = coord(:,2) + 436;
    
    % perform transformation, load matlab transform object
    %(transforms were computed in micron space)
    
    load('/Users/ashwin/Google Drive/Zfish/LowEMtoHighEM/tformRough-LowEMtoHighEM-set2-Elavl3-Mnx-SB-set4.mat');
    
    % transform from voxels to microns.
    coord(:,1) = coord(:,1) * 0.798 ;
    coord(:,2) = coord(:,2) * 0.798 ;
    coord(:,3) = coord(:,3) * 2 ;
    
    % transform from NG space to z-brainAtlas space using a precomputed
    % transfromation matrix that was loaded above.
    
    pointCloudZB = transformPointsForward(tformRough, coord);
    
    
elseif mipLevel ==4
    
    coord = pointCloudNG ;
    
    % compute offset as data in NG is offset from origin
    offset = [920,752,16400] .* [80, 80, 45]; % offset at mip4 number can be obtained from NG .info file
    
    
    coord(:,1) = coord(:,1) - offset(1);
    coord(:,2) = coord(:,2) - offset(2);
    coord(:,3) = coord(:,3) - offset(3);
    
    % adjust voxel size to LM space
    % voxel --> micron
    
    coord = coord ./ 1000;
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
    
    % perform transformation load matlab transform object to being to
    % low-res space
    load('tformRough-LowEMtoHighEM-set2.mat')
    
    % transform from voxels to microns
    coord(:,1) = coord(:,1) * 0.798 ;
    coord(:,2) = coord(:,2) * 0.798 ;
    coord(:,3) = coord(:,3) * 2 ;
    
    % transform from NG space to Atlas space using a precomputed
    % transfromation matrix.
    
    pointCloudZB = transformPointsForward(tformRough, coord);
end
end



