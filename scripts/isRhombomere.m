function [inBoundary] = isRhombomere(cellID,somaLocations)
% isRhombomere deteremines which rhombomere the root node of cellID
% lies in.
% inBoundary is nx6 where (n,1) is the cellID
% (n,2:6) are r3:r7


if ismac
    addpath(genpath('/Users/ashwin/Documents/LabWork'));
    fname  = '/Users/ashwin/Google Drive/Zfish/LowEMtoHighEM//SWC_all/combinedConsensus-resampled/';
    load('/Users/ashwin/Google Drive/Zfish/RefBrains/MaskDatabase.mat');
    imagesFile = '/Users/ashwin/Google Drive/Zfish/RefBrains/SpinalBackfills.tif';
    
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    fname =  '/usr/people/ashwinv/seungmount/research/Ashwin/Z-Brain/LowEMtoHighEM/SWC_all/consensus-20180920/swc/'
    load('/usr/people/ashwinv/seungmount/research/Ashwin/Z-Brain/RefBrains/MaskDatabase.mat');
    imagesFile = '/usr/people/ashwinv/seungmount/research/Ashwin/Z-Brain/RefBrains/SpinalBackfills.tif';
end


% need these params to set the global scale
imagesInfo = imfinfo(imagesFile);
width = imagesInfo(1).Width;
height = imagesInfo(1).Height;
imagesRef = imref3d([width,height,length(imagesInfo)],0.798,0.798,2); % set the correct dimensions for the images using the ref object


IdsHighlight = 221:225;        % r3:r7 ;

for i = 1:length(IdsHighlight)
    maskImageHighlight(i).Ids = IdsHighlight(i);
    maskImageHighlight(i).name = MaskDatabaseNames(IdsHighlight(i));
    maskImageHighlight(i).image = reshape(full(MaskDatabase(:,IdsHighlight(i))),height,width,length(imagesInfo));
end

%  load and transfrom skeleton to z-brain atlas space
if nargin == 1
    
    rootNodePlane = zeros(size(cellID,1),3);
    %swc_zbrian = zeros(size(cellID,1),1);
    
    if length(cellID) == 1
        swc_zbrian = SwctoZbrian(cellID);
        rootNodePlane(1,:) = [swc_zbrian{1}.X(1),swc_zbrian{1}.Y(1),swc_zbrian{1}.Z(1)]; % get root node coords in Z-brain space.
    else
        swc_zbrian = SwctoZbrian(cellID);
    end
    
    for i = 1:length(cellID)
        if ~isempty(swc_zbrian{i})
            rootNodePlane(i,:) = [swc_zbrian{i}.X(1),swc_zbrian{i}.Y(1),swc_zbrian{i}.Z(1)]; % get root node coords in Z-brain space.
        else
            rootNodePlane(i,:) = NaN;
        end
    end
    
    
    % determine which rhombomere the rootnode lies in.
    
    inBoundary = struct('cellID',zeros(length(cellID),1),'rootNode',zeros(length(cellID),3),'r3',zeros(length(cellID),1),'r4',zeros(length(cellID),1), ...
        'r5', zeros(length(cellID),1), 'r6', zeros(length(cellID),1), 'r7',zeros(length(cellID),1),'none',zeros(length(cellID),3));
    
    for j = 1:length(cellID)
        if ~isnan(rootNodePlane(j,:))
            inBoundary.cellID(j) = cellID(j);
            inBoundary.rootNode(j,:) = rootNodePlane(j,:);
            for i = 1:length(IdsHighlight)
                invertedPlaneinMicrons = imagesRef.ImageExtentInWorldZ-rootNodePlane(j,3);
                invertedPlaneinVoxel = round(invertedPlaneinMicrons/imagesRef.PixelExtentInWorldZ);
                [B,L] = bwboundaries(maskImageHighlight(i).image(:,:,invertedPlaneinVoxel));
                if isempty(B) ~=1
                    inveratedPlaneinMicrons = rootNodePlane(j,3);
                    temp = cell2mat(B);
                    boundaries(:,1) =  temp(:,1).*[imagesRef.PixelExtentInWorldX];
                    boundaries(:,2) = temp(:,2).*[imagesRef.PixelExtentInWorldY];
                    boundaries(:,3) = invertedPlaneinMicrons*ones(size(temp,1),1);
                    if i ==1
                        inBoundary.r3(j) = inpolygon(rootNodePlane(j,1),rootNodePlane(j,2),boundaries(:,2), boundaries(:,1));
                        clear boundaries;
                    elseif i ==2
                        inBoundary.r4(j) = inpolygon(rootNodePlane(j,1),rootNodePlane(j,2),boundaries(:,2), boundaries(:,1));
                        clear boundaries;
                    elseif i ==3
                        inBoundary.r5(j) = inpolygon(rootNodePlane(j,1),rootNodePlane(j,2),boundaries(:,2), boundaries(:,1));
                        clear boundaries;
                    elseif i ==4
                        inBoundary.r6(j) = inpolygon(rootNodePlane(j,1),rootNodePlane(j,2),boundaries(:,2), boundaries(:,1));
                        clear boundaries;
                    elseif i ==5
                        inBoundary.r7(j) = inpolygon(rootNodePlane(j,1),rootNodePlane(j,2),boundaries(:,2), boundaries(:,1));
                        clear boundaries;
                    end
                else
                    inBoundary.none(j,:) = rootNodePlane(j,:);
                end
            end
        else
            inBoundary.cellID(j) = cellID(j);
        end
    end
    
else
    
    % determine which rhombomere the rootnode lies in.
    rootNodePlane = somaLocations;
    numCells = size(somaLocations,1);
    
    inBoundary = struct('rootNode',zeros(numCells,3),'r3',zeros(numCells,1),'r4',zeros(numCells,1), ...
        'r5', zeros(numCells,1), 'r6', zeros(numCells,1), 'r7',zeros(numCells,1),'none',zeros(numCells,3));
    
    for j = 1:numCells
        if ~isnan(rootNodePlane(j,:))
            %inBoundary.cellID(j) = cellID(j);
            inBoundary.rootNode(j,:) = rootNodePlane(j,:);
            for i = 1:length(IdsHighlight)
                invertedPlaneinMicrons = rootNodePlane(j,3);
                invertedPlaneinVoxel = imagesRef.ImageSize(3)-round(invertedPlaneinMicrons/imagesRef.PixelExtentInWorldZ);
                [B,L] = bwboundaries(maskImageHighlight(i).image(:,:,invertedPlaneinVoxel));
                if isempty(B) ~=1
                    temp = cell2mat(B);
                    boundaries(:,1) =  temp(:,1).*[imagesRef.PixelExtentInWorldX];
                    boundaries(:,2) = temp(:,2).*[imagesRef.PixelExtentInWorldY];
                    boundaries(:,3) = invertedPlaneinMicrons*ones(size(temp,1),1);
                    if i ==1
                        inBoundary.r3(j) = inpolygon(rootNodePlane(j,1),rootNodePlane(j,2),boundaries(:,2), boundaries(:,1));
                        clear boundaries;
                    elseif i ==2
                        inBoundary.r4(j) = inpolygon(rootNodePlane(j,1),rootNodePlane(j,2),boundaries(:,2), boundaries(:,1));
                        clear boundaries;
                    elseif i ==3
                        inBoundary.r5(j) = inpolygon(rootNodePlane(j,1),rootNodePlane(j,2),boundaries(:,2), boundaries(:,1));
                        clear boundaries;
                    elseif i ==4
                        inBoundary.r6(j) = inpolygon(rootNodePlane(j,1),rootNodePlane(j,2),boundaries(:,2), boundaries(:,1));
                        clear boundaries;
                    elseif i ==5
                        inBoundary.r7(j) = inpolygon(rootNodePlane(j,1),rootNodePlane(j,2),boundaries(:,2), boundaries(:,1));
                        clear boundaries;
                    end
                    
                else
                    inBoundary.none(j,:) = rootNodePlane(j,:);
                end
            end
        else
            inBoundary.cellID(j) = cellID(j);
        end
    end
end
end

