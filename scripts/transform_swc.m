function [imgXZ, imgYZ] =  transform_swc_AV(cellID, neuronColor, IdsHighlight, displayBoundary, displayRefBrain, savePref)
% transfroms cellID fronm NG space to Z-brain space.
% cellID is a nx1 vector
% neuroncolor is [r,g,b] color of the neuron
% IdsHilight are the brain regions from the maskdatabase that need to be
% hilighted
% displayRefBrain is true if the background needs to be an image plane from
% the zbrain atlas.

if nargin<6
    savePref = [];
end

if isempty(cellID)
    disp('********no cells here**********');
    return;
end

if size(neuronColor,1)>1
    somataColor = neuronColor;
else
    neuronColor = repmat(neuronColor,(numel(cellID)),1);
    somataColor = neuronColor;
end

% add path .swc files of reconstructed neurons

fname  = '../skeletons/combinedConsensus-resampled/';


%% Read reference brain or any brain from the ref atlas
% NOTE: Original files needs to be rotated ccw 90 and the zaxis needs to be
% reversed in FIJI only.

%imageFileName = 'aLZBTm90-37-B45-60-footPrint-rightPos.tif';
%imageFileName = 'aLZBPointsTm90-37-B45-60-norm-COM-RightPos.tif'
%imageFileName = 'RightPos_2020_02_27_v1.tif';
%imageFileName = 'ZBB_gad1b-GFP.tif'
imageFileName = 'ZBB_ml2_mnx1-GFP.tif';
%imageFileName = 'mbrainTransformedWarp-Purk.tif';

if ismac
    imageFilePath = '/Users/ashwin/Google Drive/Zfish/RefBrains/';
else
    imageFilePath = '/usr/people/ashwinv/seungmount/research/Ashwin/Z-Brain/RefBrains/';
end

imagesInfo = imfinfo(fullfile(imageFilePath,imageFileName));
width = imagesInfo(1).Width;
height = imagesInfo(1).Height;
%images = zeros(height, width,length(imagesInfo),'uint8');
images = zeros(height, width,length(imagesInfo));

if displayRefBrain
    TifLink = Tiff(fullfile(imageFilePath,imageFileName),'r');
    
    for i = 1:length(imagesInfo)
        TifLink.setDirectory(i);
        images(:,:,i) = TifLink.read();
    end
    TifLink.close();
end
%warning('off','last');
imagesRef = imref3d(size(images),0.798,0.798,2); % set the correct dimensions for the images using the ref object

%% Mask Image
% this requires a local copy of the MaskDatabase.mat file from the Zbrain atlas
% https://zebrafishexplorer.zib.de/download

if displayBoundary
    
    disp('Load Mask database'); % mask database can be obtained from the Z-Brain atlas.
    
    if ismac
        load('/Users/ashwin/Google Drive/ZFish/RefBrains/MaskDatabase.mat'); % path to the MaskDatabase.mat
    else
        load('/usr/people/ashwinv/seungmount/research/Ashwin/Z-Brain/RefBrains/MaskDatabase.mat');
    end
    
    IdsBase = [114,221:225]; %  hardcoded to include the mask for rhombomeres r3-7
    
    for i = 1:length(IdsBase)
        maskImage(i).Ids = IdsBase(i);
        maskImage(i).name = MaskDatabaseNames(IdsBase(i));
        maskImage(i).image = reshape(full(MaskDatabase(:,IdsBase(i))),height,width,length(imagesInfo));
    end
end

% other masks that need to be displayed!

if IdsHighlight
    
    %IdsHighlight = [238,235,186];        % Vestibular Clusters
    %Ids =  [Ids,184,185,187:190];        % Reticulospinal Clusters
    %IdsHighlight =  [135,186,246,250];   % Gad1b-s2, Gly2-s2, VglutS3, vglut2-strip4
    %Ids =  [93,131,130,134];             % Cerebellum
    IdsHighlight = [91];
    
    for i = 1:length(IdsHighlight)
        maskImageHighlight(i).Ids = IdsHighlight(i);
        maskImageHighlight(i).name = MaskDatabaseNames(IdsHighlight(i));
        maskImageHighlight(i).image = reshape(full(MaskDatabase(:,IdsHighlight(i))),height,width,length(imagesInfo));
    end
end


%% get swc coordinate

cellID = cellID(isExistReRoot(cellID));
swc = cell(1,size(cellID,1));
swc_new = cell(1,size(cellID,1));
RCdist = [];
MLdist = [];
DVdist = [];

for i = 1:length(cellID)
    filename = sprintf('%d_reRoot_reSample_5000.swc',cellID(i));
    if exist(fullfile(fname,filename))
        fID = fopen(fullfile(fname,filename));
        swc{i} = textscan(fID, '%f %f %f %f %f %f %f','HeaderLines',6,'CollectOutput',true);
        fclose(fID);
    elseif exist(fullfile(fname,sprintf('%d.swc',cellID(i))))
        fID = fopen(fullfile(fname,sprintf('%d.swc',cellID(i))));
        swc{i} = textscan(fID, '%f %f %f %f %f %f %f','HeaderLines',0,'CollectOutput',true);
        fclose(fID);
    else
        continue;
    end
    clear coord;
    coord = swc{i}{1,1}(:, 3:5) ;
    
    % compute offset as data in NG is offset from origin
    
    offset = [920,752,16400] .* [80, 80, 45]; % offset at mip4 number can be obtained from NG .info file
    
    coord(:,1) = coord(:,1) - offset(1);
    coord(:,2) = coord(:,2) - offset(2);
    coord(:,3) = coord(:,3) - offset(3);
    
    % adjust voxel size to LM space
    %voxel --> micron
    
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
    
    % perform transformation
    % load the Matlba transform object to transform from NG space to Zbrain
    % space
    

    load('../matFiles/tformRough-LowEMtoHighEM-set2-Elavl3-Mnx-SB-set4.mat');
  
    
    % transform from voxels to microns
    
    coord(:,1) = coord(:,1) * 0.798 ;
    coord(:,2) = coord(:,2) * 0.798 ;
    coord(:,3) = coord(:,3) * 2 ;
    
    % transform from NG space to Atlas space using a precomputed
    % transfromation matrix.
    
    coord_transformed = transformPointsForward(tformRough, coord);
    swc_new{i}(:,3:5) = coord_transformed;
    % dlmwrite('76748_LM.swc',  swc_new, ' ');
    
    % compute edges
    RCdist = [RCdist;swc_new{1,i}(:,4)];
    MLdist = [MLdist;swc_new{1,i}(:,3)];
    DVdist = [DVdist;swc_new{1,i}(:,5)];
    
    % get mask plane and image (in mircron space, atlas space) through which the root node traverses.
    
    rootNodePlaneXY(i) = round(coord_transformed(1,3)); % in micron space
    rootNodePlaneXZ(i) = round(coord_transformed(1,2));
    rootNodePlaneYZ(i) = round(coord_transformed(1,1));
    
    
    %[~, rootNodePlaneindex(i)] = min(coord_transformed(:,1)); % in micron space
    % rootNodePlane(i) = round(coord_transformed(rootNodePlaneindex(i),3));
    if displayRefBrain
        rootNodeImageXY(:,:,i) = images(:,:,round(rootNodePlaneXY(i)/imagesRef.PixelExtentInWorldZ));
    end
end
%% construct a grided volume

% deteremine the mean plane through which the root nodes traverses

meanRootNodePlaneXY = median(rootNodePlaneXY);
meanRootNodePlaneXZ = median(rootNodePlaneXZ);
meanRootNodePlaneYZ = median(rootNodePlaneYZ);


% sample above and below root node plane

%meanRootNodePlaneXY = mean(rootNodePlaneXY(rootNodePlaneXY>quantile(rootNodePlaneXY,0.5) & rootNodePlaneXY<quantile(rootNodePlaneXY,0.95)));
%meanRootNodePlaneXY = meanRootNodePlaneXY;


% determine the mean plane image
meanImageXY = images(:,:,round(meanRootNodePlaneXY/imagesRef.PixelExtentInWorldZ));

for i = 1:size(images,3)
    meanImageYZ(i,:)  = images(:,round(meanRootNodePlaneYZ/imagesRef.PixelExtentInWorldX),i);
end



if displayRefBrain
    
    for ii = 1:size(rootNodePlaneXZ,2)
        for i = 1:size(images,3)
            rootNodeImageXZ(i,:,ii)  = images(round(rootNodePlaneXZ(ii)/imagesRef.PixelExtentInWorldX),:,i);
        end
    end
    
    
    for ii = 1:size(rootNodePlaneYZ,2)
        for i = 1:size(images,3)
            rootNodeImagYZ(i,:,ii)  = images(:,round(rootNodePlaneYZ(ii)/imagesRef.PixelExtentInWorldX),i);
        end
    end
    
    
    
    % mean projection
    meanImageXY = mean(rootNodeImageXY,3);
    meanImageXZ = mean(rootNodeImageXZ,3);
    meanImageYZ = mean(rootNodeImagYZ,3);
    %
    
    % Maximum intensity projections
    %    meanImageXY = max(images(:,:,round(min(rootNodePlaneXY)/imagesRef.PixelExtentInWorldZ):round(max(rootNodePlaneXY)/imagesRef.PixelExtentInWorldZ)),[],3);
    %meanImageXY = max(images(:,:,:),[],3);
    
    %     % max projection
    %         meanImageXY = max(rootNodeImageXY,[],3);
    %         meanImageXZ = max(rootNodeImageXZ,[],3);
    %         meanImageYZ = max(rootNodeImagYZ,[],3);
    %
    % resize the mean image to micron space
    imgXY = imresize(meanImageXY,[imagesRef.ImageExtentInWorldY,imagesRef.ImageExtentInWorldX]);
    imgXZ = imresize(meanImageXZ,[imagesRef.ImageExtentInWorldZ,imagesRef.ImageExtentInWorldX]);
    
    %imgXZ = imcomplement(medfilt2(imadjust(imgXZ/max(imgXZ(:)))));
    imgXZ = imcomplement(medfilt2(imgXZ));
    
    
    imgYZ = imresize(meanImageYZ,[imagesRef.ImageExtentInWorldZ,imagesRef.ImageExtentInWorldY]);
    %imgYZ = imcomplement(medfilt2(imadjust(imgYZ/max(imgYZ(:)))));
    imgYZ = imcomplement(medfilt2(imgYZ));
end

%imgXY = flip(imgXY,2);

% construct the mean plane as a surface

if displayRefBrain
    
    
    figure('Units','normalized','Position',[0,0,1,1]);
    % XY
    [X,Y] = meshgrid(linspace(0,ceil(imagesRef.ImageExtentInWorldX), ceil(imagesRef.ImageExtentInWorldX)),...
        linspace(0,ceil(imagesRef.ImageExtentInWorldY),ceil(imagesRef.ImageExtentInWorldY)));
    Zplane = 276*ones(size(X)); % this is to render the image layer to the back for visualization
    hsurf1 = surface(X,Y,Zplane,imcomplement(medfilt2(imadjust(imgXY/max(imgXY(:))))),'FaceColor','flat','EdgeColor','none');
    colormap(gray);
    %hsurf1 = surface(X,Y,Zplane,abs(medfilt2(imgXY)),'FaceColor','flat','EdgeColor','none');
    %hsurf1 = surface(X,Y,Zplane,imcomplement(imgXY),'FaceColor','flat','EdgeColor','none');
    
    %colormap(cmapRed);
    %
    lighting gouraud;
    material shiny;
    %alpha(hsurf1, 0.8);
    clear X;
    clear Y;
    
end

daspect([1,1,1]);
hold on;

% IDs to be highlighted (e.g. Tangential Vestibular Nucleus)

if IdsHighlight
    index = 0.1;
    
    for i = 1:length(IdsHighlight)
        invertedPlaneinMicrons = imagesRef.ImageExtentInWorldZ-meanRootNodePlaneXY;
        invertedPlaneinVoxel = round(invertedPlaneinMicrons/imagesRef.PixelExtentInWorldZ);
        [B,L] = bwboundaries(maskImageHighlight(i).image(:,:,invertedPlaneinVoxel));
        %invertedPlaneinMicrons = meanRootNodePlane;
        %boundarySize = size(B,1);
        if size(B,1)>0
            boundaries = B{1};
            boundaries = boundaries.*[imagesRef.PixelExtentInWorldX,imagesRef.PixelExtentInWorldY];
            boundaries(:,3) = meanRootNodePlaneXY*ones(size(boundaries,1),1);
            %fill3(boundaries(:,2),boundaries(:,1),boundaries(:,3),cols(i,:),'FaceAlpha',0.3,'LineStyle','none','LineWidth',0.5);
            plot3(boundaries(:,2),boundaries(:,1),boundaries(:,3),'Color',cols(i,:),'LineWidth',1,'LineStyle','--');
            clear boundaries;
        end
        if size(B,1)>1
            boundaries = B{2};
            boundaries = boundaries.*[imagesRef.PixelExtentInWorldX,imagesRef.PixelExtentInWorldY];
            boundaries(:,3) = meanRootNodePlaneXY*ones(size(boundaries,1),1);
            %fill3(boundaries(:,2),boundaries(:,1),boundaries(:,3),cols(i,:),'FaceAlpha',0.3,'LineStyle','none','LineWidth',0.5);
            plot3(boundaries(:,2),boundaries(:,1),boundaries(:,3),'Color',cols(i,:),'LineWidth',1,'LineStyle','--');
            annotation('textbox',[0.42,0.4-index,0.5,0],'EdgeColor','none','String',maskImageHighlight(i).name,'Color',cols(i,:));
        end
        index = index+0.01;
    end
end

% plot cells in cellID

for i = 1:length(cellID)
    filename = sprintf('%d_reRoot_reSample_5000.swc',cellID(i));
    if exist(fullfile(fname,filename))
        tree = load_tree(fullfile(fname,filename));
        % if neurites are not needed comment the below three lines
        [I,J] = ind2sub(size(tree.dA),find(tree.dA));
        line([swc_new{i}(J,3) swc_new{i}(I,3)]',[swc_new{i}(J,4) swc_new{i}(I,4)]',[swc_new{i}(J,5) swc_new{i}(I,5)]',...
            'Color',[neuronColor(i,:),0.5],'LineWidth',1);
        hold on;
        % handle for soma display properties
        h = scatter3(swc_new{i}(1,3), swc_new{i}(1,4), swc_new{i}(1,5),100,'o','MarkerEdgeColor','k',...
            'MarkerFaceColor','none','LineWidth',0.2);
        hMarker = h.MarkerHandle;
        hMarker.FaceColorType = 'truecoloralpha';
        hMarker.FaceColorData = uint8(255*[somataColor(i,1);somataColor(i,2);somataColor(i,3);0.3]);

        clear I;
        clear J;
        clear tree;
    elseif exist(fullfile(fname,sprintf('%d.swc',cellID(i))))
        tree = load_tree(fullfile(fname,sprintf('%d.swc',cellID(i))));
        [I,J] = ind2sub(size(tree.dA),find(tree.dA));
        line([swc_new{i}(J,3) swc_new{i}(I,3)]',[swc_new{i}(J,4) swc_new{i}(I,4)]',[swc_new{i}(J,5) swc_new{i}(I,5)]',...
            'Color',[neuronColor(i,:),0.1],'LineWidth',0.1);
        hold on;
        h = scatter3(swc_new{i}(1,3), swc_new{i}(1,4), swc_new{i}(1,5),50,'MarkerFaceColor',somataColor(i,:),'MarkerEdgeColor','k','LineWidth',1.0);
        hMarker = h.MarkerHandle;
        hMarker.FaceColorData = uint8(255*[somataColor(i,1);somataColor(i,2);somataColor(i,3)]);
        
        clear I;
        clear J;
        clear tree;
    else
        continue;
    end
end
h = gca;

% draw boundary of the anatomical region of interst

%All masks are in voxel dimensions, need to covert them to micron space.
%masks are also inverted, meaning from V-->D.

% Base level maske IDs (used for oveall orientation, e.g. rhombomere
% boundaries)

if displayBoundary
    cols = cbrewer('qual','Accent',length(IdsBase));
    
    for i = 1:length(IdsBase)
        invertedPlaneinMicrons = imagesRef.ImageExtentInWorldZ-meanRootNodePlaneXY;
        invertedPlaneinVoxel = round(invertedPlaneinMicrons/imagesRef.PixelExtentInWorldZ);
        [B,L] = bwboundaries(maskImage(i).image(:,:,invertedPlaneinVoxel));
        %invertedPlaneinMicrons = meanRootNodePlane;
        %boundarySize = size(B,1);
        if size(B,1)>0
            boundaries = B{1};
            boundaries = boundaries.*[imagesRef.PixelExtentInWorldX,imagesRef.PixelExtentInWorldY];
            boundaries(:,3) = meanRootNodePlaneXY*ones(size(boundaries,1),1);
            %fill3(boundaries(:,2),boundaries(:,1),boundaries(:,3),cols(i,:),'FaceAlpha',0.3,'LineStyle','none','LineWidth',0.5);
            plot3(boundaries(:,2),boundaries(:,1),boundaries(:,3),'Color',[0.5,0.5,0.5],'LineWidth',2,'LineStyle','-');
            clear boundaries;
        end
        if size(B,1)>1
            boundaries = B{2};
            boundaries = boundaries.*[imagesRef.PixelExtentInWorldX,imagesRef.PixelExtentInWorldY];
            boundaries(:,3) = meanRootNodePlaneXY*ones(size(boundaries,1),1);
            %fill3(boundaries(:,2),boundaries(:,1),boundaries(:,3),cols(i,:),'FaceAlpha',0.3,'LineStyle','none','LineWidth',0.5);
            plot3(boundaries(:,2),boundaries(:,1),boundaries(:,3),'Color',[0.5,0.5,0.5],'LineWidth',2,'LineStyle',':');
        end
    end
end


set(gca, 'BoxStyle','full','YDir','reverse','ZDir','reverse','color','none');
%set(gca, 'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
axis on;
%set(gca,'ZLim',[0,imagesRef.ImageExtentInWorldZ], 'XLim',[0,imagesRef.ImageExtentInWorldX], 'YLim', [0,imagesRef.ImageExtentInWorldY]);
set(gca,'ZLim',[0,imagesRef.ImageExtentInWorldZ], 'XLim',[250,imagesRef.ImageExtentInWorldX], 'YLim', [400,800]);

%draw scale bar 1px = 0.798 us
%line([400,463],[750,750],'Color','k','LineWidth',4);
%text(500,1300,'50um','FontName','Arial','FontSize',10);
%axis off;

if ~isempty(savePref)
    folderName = sprintf('/Users/ashwin/Desktop/%s',savePref);
    mkdir(folderName);
    
    % XY view
    %set(gcf, 'Renderer','painters');
    fileName = sprintf('%s/%s-XY.png',folderName,savePref);
    export_fig(fileName,'-r300','-transparent');
    
    pause(2);
    
    % XZ
    %set(gcf, 'Renderer','painters');
    fileName = sprintf('%s/%s-XZ.png',folderName,savePref);
    view([0,0]);
    export_fig(fileName,'-r300','-transparent');
    
    %YX
    %set(gcf, 'Renderer','painters');
    fileName = sprintf('%s/%s-YZ.png',folderName,savePref);
    view([90,0]);
    export_fig(fileName,'-r300','-transparent');
    
    close all
    
    figure('Units','normalized','Position',[0,0,1,1]);
    imagesc(imcomplement(abs(imgXZ)));
    colormap(gray);
    daspect([1,1,1]);
    fileName = sprintf('%s/%s-BackgroundXZ.png',folderName,savePref);
    export_fig(fileName,'-r300','-transparent');
    close all
    
    figure('Units','normalized','Position',[0,0,1,1]);
    imagesc(imcomplement(abs(imgYZ)));
    colormap(gray);
    daspect([1,1,1]);
    fileName = sprintf('%s/%s-BackgroundYZ.png',folderName,savePref);
    export_fig(fileName,'-r300','-transparent');
    close all;
    
end
%title(sprintf('Zbr plane: %1d',138-round(meanRootNodePlaneXY/2)));
%title(sprintf('cellID: %1d',cellID));
%title(sprintf('Zbr dorsal plane: %1d \n Zbr ventral plane: %1d', round(min(rootNodePlane(rootNodePlane~=0))/2), round(max(rootNodePlane)/2)));

end

