% plots for clustering data

load ('../matFiles/MatOrder_CO_top500_2blocks_gamma038_08062020.mat');
load ('../matFiles/cellIDType_CO_top500_2blocks_gamma038_08062020.mat');

cellIDType = cellstr(cellIDType_CO_top500_2blocks_gamma038_08062020);

block1_cleaned.cellIDs = MatOrder_CO_top500_2blocks_gamma038_08062020(strcmp(cellIDType,'_Axl_'));
block2_cleaned.cellIDs = [MatOrder_CO_top500_2blocks_gamma038_08062020(strcmp(cellIDType,'_Int_')),...
                          MatOrder_CO_top500_2blocks_gamma038_08062020(strcmp(cellIDType,'_DOs_'))];

block1_cleaned.color = [1,0.5,0];
block2_cleaned.color = [0,0.5,1];

load ('../matFiles/MatOrder_CO_OM_top500_2blocks_gamma038_08072020.mat');
load ('../matFiles/cellIDType_CO_OM_top500_2blocks_gamma038_08072020.mat');

cellIDType_OM = cellstr(cellIDType_CO_OM_top500_2blocks_gamma038_08072020);
block1_cleaned_submod.cellIDs = MatOrder_CO_OM_top500_2blocks_gamma038_08072020(strcmp(cellIDType_OM,'mod_1'));
block2_cleaned_submod.cellIDs = MatOrder_CO_OM_top500_2blocks_gamma038_08072020(strcmp(cellIDType_OM,'mod_2'));


temp = distinguishable_colors(10,'w');
block1_cleaned_submod.color = temp(9,:);
block2_cleaned_submod.color = temp(8,:);

load ('../matFiles/AllCells.mat');
load ('../matFiles/ConnMatrixPre_cleaned.mat');


%LoadDataFrame
load ('../matFiles/df_cleaned.mat');

confirmedALX = [76181 76201 76187 76184 76192 76197];
confirmedDBX = [76199 76200 76182 76183 76185 76186 76189 76191 76188 ];
confirmedBARHL = [76198 76190 76193 76194 76195 76196 ];
confirmedIntegrators = [confirmedALX,confirmedDBX,confirmedBARHL];


% ABD motor neurons
ABDr_CellIDs = [82140,82145,82143,77648,82146,77302,77301,82194,77705,77710,77305,77672,77300,77709,77661];
ABDc_CellIDs = [82213,77654,77646,82212,77295,81172,77657,82197,77682,77154, 77652,77658,77628,77292,77688,82195,77296];
ABDIr_CellIDs = [78553,77150,77886,78547,77631,77158,77665,78552,77668,77618,77634,78556];
ABDIc_CellIDs = [78574,79051,77148,77625,77144,77692,77641,79066,77640,77643];


%%

block2_cleaned_submod.noVest = block2_cleaned_submod.cellIDs(~isVestibular(block2_cleaned_submod.cellIDs));
block1_cleaned_submod.noVest = block1_cleaned_submod.cellIDs(~isVestibular(block1_cleaned_submod.cellIDs));

OMcells = [block2_cleaned_submod.cellIDs(~isVestibular(block2_cleaned_submod.cellIDs)),...
    block1_cleaned_submod.cellIDs(~isVestibular(block1_cleaned_submod.cellIDs)),...
    block1_cleaned_submod.cellIDs(isVestibular(block1_cleaned_submod.cellIDs)),...
    block2_cleaned_submod.cellIDs(isVestibular(block2_cleaned_submod.cellIDs))];

ABD_cells =    [ABDr_CellIDs,ABDc_CellIDs,ABDIr_CellIDs,ABDIc_CellIDs];

[~,OMcellsIndex] = ismember(OMcells,AllCells);
[~,ABDcellsIndex] = ismember(ABD_cells,AllCells);

gap = 20;

tempMat = [];
matSize = size(OMcells,2)+gap+size(ABD_cells,2);
tempMat = zeros(matSize);

tempMat(1:size(OMcells,2),1:size(OMcells,2)) = ConnMatrixPre_cleaned(OMcellsIndex,OMcellsIndex);
tempMat(size(OMcells,2)+gap+1:end,1:size(OMcells,2)) = ConnMatrixPre_cleaned(ABDcellsIndex,OMcellsIndex);

tempMatThreshold = [];
tempMatThreshold = tempMat;
tempMatThreshold(tempMatThreshold>=5) =5;

figure;
cspy(tempMatThreshold,'Colormap',colorcet('R3','N',5),'Levels',5,'MarkerSize',8);
c = colorbar;
c.Ticks = [1:5];
c.TickLabels = [1:5];
hold on;

line([0,size(tempMat,2)],[size(block2_cleaned_submod.noVest,2),size(block2_cleaned_submod.noVest,2)],'color','k');
line([size(block2_cleaned_submod.noVest,2),size(block2_cleaned_submod.noVest,2)],[0,size(tempMat,2)],'color','k');
line([0,size(tempMat,2)],[size(block2_cleaned_submod.noVest,2),size(block2_cleaned_submod.noVest,2)]+ size(block1_cleaned_submod.noVest,2),'color','k');
line([size(block2_cleaned_submod.noVest,2),size(block2_cleaned_submod.noVest,2)]+size(block1_cleaned_submod.noVest,2),[0,size(tempMat,2)],'color','k');

line([0,size(OMcells,2)],[size(OMcells,2),size(OMcells,2)],'color','k');
line([size(OMcells,2),size(OMcells,2)],[0,size(OMcells,2)],'color','k');
line([0,size(OMcells,2)],[size(OMcells,2),size(OMcells,2)]+gap,'color','k');
line([0,size(OMcells,2)],[size(OMcells,2),size(OMcells,2)]+gap+32,'color','k');

box on;

%% Pre and Post locations (long time)

% block1

block1_cleaned.cellIDs(block1_cleaned.cellIDs ==76202) = []; % remove Mauthner for further analysis.
block1_cleaned.dendSites = [];
block1_cleaned.dendpartnerSites = [];
for i = 1:length(block1_cleaned.cellIDs)
    [temp,preIDs] = SynapticPartners(block1_cleaned.cellIDs(i),1,df_cleaned);
    preIDs = preIDs(temp~=block1_cleaned.cellIDs(i)); % remove self touches
    block1_cleaned.dendSites = [block1_cleaned.dendSites;PrePartnerCoordinates(preIDs,df_cleaned)];
    block1_cleaned.dendpartnerSites = [block1_cleaned.dendpartnerSites;PostPartnerCoordinates(preIDs,df_cleaned)];
    clear preIDs;
    clear temp;
    i
end
block1_cleaned.dendSites = TransformPoints(block1_cleaned.dendSites,0);
block1_cleaned.dendpartnerSites = TransformPoints(block1_cleaned.dendpartnerSites,0);


block2_cleaned.dendSites = [];
block2_cleaned.dendpartnerSites = [];
for i = 1:length(block2_cleaned.cellIDs)
    [temp,preIDs] = SynapticPartners(block2_cleaned.cellIDs(i),1,df_cleaned);
    preIDs = preIDs(temp~=block2_cleaned.cellIDs(i)); % remove self touches
    block2_cleaned.dendSites = [block2_cleaned.dendSites;PrePartnerCoordinates(preIDs,df_cleaned)];
    block2_cleaned.dendpartnerSites = [block2_cleaned.dendpartnerSites;PostPartnerCoordinates(preIDs,df_cleaned)];
    clear preIDs;
    i
end
block2_cleaned.dendSites = TransformPoints(block2_cleaned.dendSites,0);
block2_cleaned.dendpartnerSites = TransformPoints(block2_cleaned.dendpartnerSites,0);

% plot somata

block1_cleaned.removeSoma = [82219,78151,78243,81524,81781,79156,80834,78926,79762,77240];
block2_cleaned.removeSoma = [78617,81066,78560,81329,80990,81021];

block1_cleaned.cleanSoma = block1_cleaned.cellIDs(~ismember(block1_cleaned.cellIDs,block1_cleaned.removeSoma));
block2_cleaned.cleanSoma = block2_cleaned.cellIDs(~ismember(block2_cleaned.cellIDs,block2_cleaned.removeSoma));

block1_cleaned.noSoma = [79407,79732,79921,79403,79559,77799,77039,80842,81519,79542,76950,77057,77084,79452,81661,76964,76931,81540,77101,78062,79280,77024,78119,78870,78847,79951,81151,78104,78143,78167,77103,77079,79687,76626,79756,78881,81147,76918,79811,79695,77082,78158,81145,81845,79580,79907,79581,81557,79227,81611,78914,79301,76922,79771,76953,79503,78121,79923,76938,79240,77054,76924,76936,78125,78123,78939,79341,77022,79734,78116,77056,77098,78877,79383,78128,78874,79560,81410,77048,78925];
block2_cleaned.noSoma = [78583,78883,77122,78406,79746,77581,76690,80548,77848,78860,77447,77872,79072,80606,77651,81374,81817,78443,78544,79055,77163,81297,77132,77433,81407,77349,80974,78351,77140,79062,77344,81537,78441,76697,81391,78629,77467,78150,78650,78226,80629,80995,76677,81431,76673,77797,79042,80539,76618,78421,79074,78297,77607,77162,79069,79022,78858,81793,80647,81317,81060,81007,79046,80757,77805,80743,80472,76701,80294,78641,79743,78649,81765,78254,77121,81322,77151,79720,79976,77329,77868,78540,80510,76622,76661,78601,81792,78543,81161,76627,81293,78246,77146,78911,78255,81774,78356,80800,81423,77434,77465,80679,77435,77341,79054,78357,77592,78949,77845,81336,81311,77591,77239,78667,78157,77124,80956,79722,77327,77844,77621,81637,81400,78646,79085,81002,77460,79080,77437,79048,77142,77389,78241,80947,78404,77636,77821,81312,79953,80681,79058,79033,77822,80746,80750,77656,80850,79067,80625,78453,76623,81363,81550,80626,81683,78633,79044,81395,77152,81443,80262,80315,80821,77602,77251,81408,77369,81580,81417,76625,78696,77580,79852,77126,81295,81559,81338];

block1_cleaned.origins = getOrigin(block1_cleaned.cleanSoma);
block2_cleaned.origins = getOrigin(block2_cleaned.cleanSoma);

block1_cleaned.inVol = block1_cleaned.origins(~ismember(block1_cleaned.cleanSoma,block1_cleaned.noSoma),:);
block2_cleaned.inVol = block2_cleaned.origins(~ismember(block2_cleaned.cleanSoma,block2_cleaned.noSoma),:);

block1_cleaned.outVol = block1_cleaned.origins(ismember(block1_cleaned.cleanSoma,block1_cleaned.noSoma),:);
block2_cleaned.outVol = block2_cleaned.origins(ismember(block2_cleaned.cleanSoma,block2_cleaned.noSoma),:);



figure;

% scatter plot for neurons with somata inside the volume.

% YZ
 h = scatterhist([block1_cleaned.origins(~ismember(block1_cleaned.cleanSoma,block1_cleaned.noSoma),2);...
     block2_cleaned.origins(~ismember(block2_cleaned.cleanSoma,block2_cleaned.noSoma),2)],...
     [block1_cleaned.origins(~ismember(block1_cleaned.cleanSoma,block1_cleaned.noSoma),3);...
     block2_cleaned.origins(~ismember(block2_cleaned.cleanSoma,block2_cleaned.noSoma),3)],...
     'Group',[ones(sum(~ismember(block1_cleaned.cleanSoma,block1_cleaned.noSoma)),1);2*ones(sum(~ismember(block2_cleaned.cleanSoma,block2_cleaned.noSoma)),1)],...
     'Kernel','on','color',[block1_cleaned.color;block2_cleaned.color],'Marker','o','MarkerSize',5,'LineWidth',2);
  set(h(1),'YDir','reverse');
  set(h(2),'Visible','on');
 set(h(3),'XDir','reverse');
 hold on;
scatter(h(1),block1_cleaned.origins(ismember(block1_cleaned.cleanSoma,block1_cleaned.noSoma),2),...
    block1_cleaned.origins(ismember(block1_cleaned.cleanSoma,block1_cleaned.noSoma),3),...
    20,'MarkerEdgeColor',block1_cleaned.color,'Marker','o','LineWidth',0.25);
scatter(h(1),block2_cleaned.origins(ismember(block2_cleaned.cleanSoma,block2_cleaned.noSoma),2),...
    block2_cleaned.origins(ismember(block2_cleaned.cleanSoma,block2_cleaned.noSoma),3),...
    20,'MarkerEdgeColor',block2_cleaned.color,'Marker','o','LineWidth',0.25);
daspect([1,1,1]);

%XY
figure;
h = scatterhist([block1_cleaned.origins(~ismember(block1_cleaned.cleanSoma,block1_cleaned.noSoma),1);...
    block2_cleaned.origins(~ismember(block2_cleaned.cleanSoma,block2_cleaned.noSoma),1)],...
    [block1_cleaned.origins(~ismember(block1_cleaned.cleanSoma,block1_cleaned.noSoma),2);...
    block2_cleaned.origins(~ismember(block2_cleaned.cleanSoma,block2_cleaned.noSoma),2)],...
    'Group',[ones(sum(~ismember(block1_cleaned.cleanSoma,block1_cleaned.noSoma)),1);2*ones(sum(~ismember(block2_cleaned.cleanSoma,block2_cleaned.noSoma)),1)],...
    'Kernel','on','color',[block1_cleaned.color;block2_cleaned.color],'Marker','o','MarkerSize',5,'LineWidth',2);
set(h(1),'YDir','reverse','YLim',[400,800],'XLim',[250,400]);
set(h(2),'Visible','on');
set(h(3),'XDir','reverse','Visible','on');
hold on;

scatter(h(1),block1_cleaned.origins(ismember(block1_cleaned.cleanSoma,block1_cleaned.noSoma),1),...
    block1_cleaned.origins(ismember(block1_cleaned.cleanSoma,block1_cleaned.noSoma),2),...
    20,'MarkerEdgeColor',block1_cleaned.color,'Marker','o','LineWidth',0.25);
scatter(h(1),block2_cleaned.origins(ismember(block2_cleaned.cleanSoma,block2_cleaned.noSoma),1),...
    block2_cleaned.origins(ismember(block2_cleaned.cleanSoma,block2_cleaned.noSoma),2),...
    20,'MarkerEdgeColor',block2_cleaned.color,'Marker','o','LineWidth',0.25);
daspect([1,1,1]);



% all somata scatter plot
figure

h = scatterhist([block1_cleaned.origins(:,1);block2_cleaned.origins(:,1)],[block1_cleaned.origins(:,2);block2_cleaned.origins(:,2)],...
     'Group',[ones(length(block1_cleaned.cleanSoma),1);2*ones(length(block2_cleaned.cleanSoma),1)],...
     'Kernel','on','color',[block1_cleaned.color;block2_cleaned.color],'Marker','o','MarkerSize',5,'LineWidth',2);
set(h(1),'YDir','reverse','YLim',[400,800],'XLim',[250,400]);
set(h(2),'Visible','on');
set(h(3),'XDir','reverse','Visible','on');
daspect([1,1,1]);


% plot histograms of the inputs
Xcords_Pre = [block1_cleaned.dendSites(1:10:end,1);block2_cleaned.dendSites(1:10:end,1)];
Ycords_Pre = [block1_cleaned.dendSites(1:10:end,2);block2_cleaned.dendSites(1:10:end,2)];
Zcords_pre = [block1_cleaned.dendSites(1:10:end,3);block2_cleaned.dendSites(1:10:end,3)];

GpIDs_Pre = [ones(length(block1_cleaned.dendSites(1:10:end,1)),1);2*ones(length(block2_cleaned.dendSites(1:10:end,1)),1)];

h = scatterhist(Xcords_Pre,Ycords_Pre,'Group',GpIDs_Pre,'Kernel','on','color',[block1_cleaned.color;block2_cleaned.color],...
    'Marker','.','MarkerSize',1,'LineWidth',2);
set(h(1),'YDir','reverse','YLim',[400,800],'XLim',[250,400]);
set(h(2),'Visible','on');
set(h(3),'XDir','reverse','Visible','on');
daspect([1,1,1]);

figure;
h = scatterhist(Ycords_Pre,Zcords_pre,'Group',GpIDs_Pre,'Kernel','on','color',[block1_cleaned.color;block2_cleaned.color],...
    'Marker','.','MarkerSize',1,'LineWidth',2);
set(h(1),'YDir','reverse');
set(h(3),'XDir','reverse');
daspect([1,1,1]);




% n =100;
% Xlin = linspace(min(block1_cleaned.dendSites(:,1)), max(block1_cleaned.dendSites(:,1)),n);
% Ylin = linspace(min(block1_cleaned.dendSites(:,2)), max(block1_cleaned.dendSites(:,2)),n);
% 
% Xr = interp1(Xlin,1:numel(Xlin),block1_cleaned.dendSites(:,1),'nearest');
% Yr = interp1(Ylin,1:numel(Ylin),block1_cleaned.dendSites(:,2),'nearest');
% 
% Z = accumarray([Xr,Yr],1,[n,n]);
% figure;
% contour(Z,'color',col1);


% h = scatterhist(Ycords_Pre,Zcords_pre,'Group',GpIDs_Pre,'Kernel','on','color',[col1;col2],...
%     'Marker','.','MarkerSize',1,'LineWidth',2);

% block2_cleaned_gamma038
block1_cleaned.axonSites = [];
block1_cleaned.axonPartnerSites = [];
for i = 1:length(block1_cleaned.cellIDs)
    [temp,preIDs] = SynapticPartners(block1_cleaned.cellIDs(i),2,df_cleaned);
    preIDs = preIDs(temp~=block1_cleaned.cellIDs(i)); % remove self touches
    block1_cleaned.axonSites = [block1_cleaned.axonSites;PostPartnerCoordinates(preIDs,df_cleaned)];
    block1_cleaned.axonPartnerSites = [block1_cleaned.axonPartnerSites;PrePartnerCoordinates(preIDs,df_cleaned)];
    clear preIDs;
    clear temp;
end
block1_cleaned.axonSites = TransformPoints(block1_cleaned.axonSites,0);
block1_cleaned.axonPartnerSites = TransformPoints(block1_cleaned.axonPartnerSites,0);

block2_cleaned.axonSites = [];
block2_cleaned.axonPartnerSites = [];
for i = 1:length(block2_cleaned.cellIDs)
    [temp,preIDs] = SynapticPartners(block2_cleaned.cellIDs(i),2,df_cleaned);
    preIDs = preIDs(temp~=block2_cleaned.cellIDs(i)); % remove self touches
    block2_cleaned.axonSites = [block2_cleaned.axonSites;PostPartnerCoordinates(preIDs,df_cleaned)];
    block2_cleaned.axonPartnerSites = [block2_cleaned.axonPartnerSites;PrePartnerCoordinates(preIDs,df_cleaned)];
    clear preIDs;
end
block2_cleaned.axonSites  = TransformPoints(block2_cleaned.axonSites ,0);
block2_cleaned.axonPartnerSites = TransformPoints(block2_cleaned.axonPartnerSites,0);

% plot histogram of the outputs
Xcords_Post = [block1_cleaned.axonSites(1:10:end,1);block2_cleaned.axonSites(1:10:end,1)];
Ycords_Post = [block1_cleaned.axonSites(1:10:end,2);block2_cleaned.axonSites(1:10:end,2)];
Zcords_Post = [block1_cleaned.axonSites(1:10:end,3);block2_cleaned.axonSites(1:10:end,3)];

figure;
 GpIDs_Post = [ones(length(block1_cleaned.axonSites(1:10:end,1)),1);2*ones(length(block2_cleaned.axonSites(1:10:end,1)),1)];
h = scatterhist(Xcords_Post,Ycords_Post,'Group',GpIDs_Post,'Kernel','on','color',[block1_cleaned.color;block2_cleaned.color],...
    'Marker','.','MarkerSize',0.5,'LineWidth',2);
set(h(1),'YDir','reverse','YLim',[400,800],'XLim',[250,400]);
set(h(2),'Visible','on');
set(h(3),'XDir','reverse','Visible','on');
daspect([1,1,1]);

figure;

h = scatterhist(Ycords_Post,Zcords_Post,'Group',GpIDs_Post,'Kernel','on','color',[block1_cleaned.color;block2_cleaned.color],...
    'Marker','.','MarkerSize',0.5,'LineWidth',2);
 set(h(1),'YDir','reverse');
 set(h(3),'XDir','reverse');
daspect([1,1,1]);


SomaSegIndex = SomasegIndex(block1_cleaned.inVol,block2_cleaned.inVol);
PostSegIndex = SomasegIndex(block1_cleaned.dendSites,block2_cleaned.dendSites);
PreSegIndex = SomasegIndex(block1_cleaned.axonSites,block2_cleaned.axonSites);

figure;
subplot(4,4,1)
heatmap(SomaSegIndex,'XDisplayLabels',{'modA','modO'},'YDisplayLabels',{'modA','modO'});
subplot(4,4,2)
heatmap(SomaSegIndex./[size(block1_cleaned.inVol,1),size(block1_cleaned.inVol,1);size(block2_cleaned.inVol,1),size(block2_cleaned.inVol,1)],'XDisplayLabels',{'modA','modO'},'YDisplayLabels',{'modA','modO'});
title('Soma segregation Index');
subplot(4,4,3)
heatmap(PostSegIndex./[size(block1_cleaned.dendSites,1),size(block1_cleaned.dendSites,1);size(block2_cleaned.dendSites,1),size(block2_cleaned.dendSites,1)],'XDisplayLabels',{'modA','modO'},'YDisplayLabels',{'modA','modO'});
title('PostSynapse segregation Index');
subplot(4,4,4)
heatmap(PreSegIndex./[size(block1_cleaned.axonSites,1),size(block1_cleaned.axonSites,1);size(block2_cleaned.axonSites,1),size(block2_cleaned.axonSites,1)],'XDisplayLabels',{'modA','modO'},'YDisplayLabels',{'modA','modO'});
title('PreSynapse segregation Index');

%% Peters Rule
% Do block1_cleaned neurons have the potential to make synapses onto the other
% block.
% block1_cleaned --> block1_cleaned (actual synapses/ potential synapses)
% block1_cleaned --> block2_cleaned_gamma038

for i = 1:length(block1_cleaned.cellIDs)
    [prePartner,prePartnerPSD] = SynapticPartners(block1_cleaned.cellIDs(i),1,df_cleaned);
    % remove self touches
    %block1_cleaned.block1_cleanedtoblock1_cleanedActual(i) = sum(ismember(prePartner,setdiff(block1_cleaned.cellIDs,block1_cleaned.cellIDs(i)))); 
    prePartnerPSD = prePartnerPSD(prePartner~=block1_cleaned.cellIDs(i));

    inputs = PrePartnerCoordinates(prePartnerPSD,df_cleaned); % location of presynaptic partner, on the axon
    inputs = TransformPoints(inputs,0);
    % distance matrix block1_cleaned -->block1_cleaned
    [distMat11] = pdist2(block1_cleaned.axonPartnerSites,inputs);
    block1_cleaned.actualSynapses11(i) = length(find(distMat11 == 0));
    block1_cleaned.potentialSynpses11_05(i) = length(find(distMat11 <= 0.5));
    block1_cleaned.potentialSynpses11_1(i) = length(find(distMat11 <= 1));
    block1_cleaned.potentialSynpses11_2(i) = length(find(distMat11 <= 2));
    block1_cleaned.potentialSynpses11_5(i) = length(find(distMat11 <= 5));
    block1_cleaned.potentialSynpses11_10(i) = length(find(distMat11 <=10));


    %distance Matrix block1_cleaned -->block2_cleaned_gamma038
    distMat12 = pdist2(block2_cleaned.axonPartnerSites,inputs);
    block1_cleaned.actualSynapses12(i) = length(find(distMat12 == 0));
    block1_cleaned.potentialSynpses12_05(i) = length(find(distMat12 <= 0.5));
    block1_cleaned.potentialSynpses12_1(i) = length(find(distMat12 <= 1));
    block1_cleaned.potentialSynpses12_2(i) = length(find(distMat12 <= 2));
    block1_cleaned.potentialSynpses12_5(i) = length(find(distMat12 <= 5));
    block1_cleaned.potentialSynpses12_10(i) = length(find(distMat12 <= 10));

    
    clear prePartner;
    clear prePartnerPSD;
    clear inputs;
    clear distMat11;
    clear distMat12;
end

% block2_cleaned_gamma038 --> block2_cleaned_gamma038 (actual synapses/ potential synapses)
% block2_cleaned_gamma038 --> block1_cleaned

for i = 1:length(block2_cleaned.cellIDs)
    [prePartner,prePartnerPSD] = SynapticPartners(block2_cleaned.cellIDs(i),1,df_cleaned);
    % remove self touches
    block2_cleaned.block2_cleaned_gamma038toblock2_cleaned_gamma038Actual(i) = sum(ismember(prePartner,setdiff(block2_cleaned.cellIDs,block2_cleaned.cellIDs(i)))); 
    prePartnerPSD = prePartnerPSD(prePartner~=block2_cleaned.cellIDs(i));
    inputs = PrePartnerCoordinates(prePartnerPSD,df_cleaned); % location of presynaptic partner, on the axon
    inputs = TransformPoints(inputs,0);
    % distance matrix block2_cleaned_gamma038 -->block2_cleaned_gamma038
    distMat22 = pdist2(block2_cleaned.axonPartnerSites,inputs);
    block2_cleaned.actualSynapses22(i) = length(find(distMat22 == 0));
    block2_cleaned.potentialSynpses22_05(i) = length(find(distMat22 <= 0.5));
    block2_cleaned.potentialSynpses22_1(i) = length(find(distMat22 <= 1));
    block2_cleaned.potentialSynpses22_2(i) = length(find(distMat22 <= 2));
    block2_cleaned.potentialSynpses22_5(i) = length(find(distMat22 <= 5));
    block2_cleaned.potentialSynpses22_10(i) = length(find(distMat22 <= 10));


    %distance Matrix block2_cleaned_gamma038 -->block1_cleaned
    distMat21 = pdist2(block1_cleaned.axonPartnerSites,inputs);
    block2_cleaned.actualSynapses21(i) = length(find(distMat21 == 0));
    block2_cleaned.potentialSynpses21_05(i) = length(find(distMat21 <= 0.5));
    block2_cleaned.potentialSynpses21_1(i) = length(find(distMat21 <= 1));
    block2_cleaned.potentialSynpses21_2(i) = length(find(distMat21 <= 2));
    block2_cleaned.potentialSynpses21_5(i) = length(find(distMat21 <= 5));
    block2_cleaned.potentialSynpses21_10(i) = length(find(distMat21 <= 10));

    
    clear prePartner;
    clear prePartnerPSD;
    clear inputs;
    clear distMat22;
    clear distMat21;
end

% figure;
% subplot(4,4,1)
% heatmap([mean(block1_cleaned.actualSynapses11),mean( block1_cleaned.actualSynapses12);mean( block2_cleaned_gamma038.actualSynapses21),mean(block2_cleaned_gamma038.actualSynapses22)],...
%     'ColorbarVisible','off');
% title('True synapses');
% 
% subplot(4,4,2)
% heatmap([mean(block1_cleaned.potentialSynpses11_2),mean(block1_cleaned.potentialSynpses12_2);mean( block2_cleaned_gamma038.potentialSynpses21_2),mean(block2_cleaned_gamma038.potentialSynpses22_2)],...
%     'ColorbarVisible','off');
% title('Potential synapses (2um)');
% 
% subplot(4,4,3)
% heatmap([mean(block1_cleaned.potentialSynpses11_5),mean(block1_cleaned.potentialSynpses12_5);mean( block2_cleaned_gamma038.potentialSynpses21_5),mean(block2_cleaned_gamma038.potentialSynpses22_5)],...
%     'ColorbarVisible','off');
% title('Potential synapses (5um)');
% 
% subplot(4,4,4)
% heatmap([mean(block1_cleaned.potentialSynpses11_10),mean(block1_cleaned.potentialSynpses12_10);mean( block2_cleaned_gamma038.potentialSynpses21_10),mean(block2_cleaned_gamma038.potentialSynpses22_10)],...
%     'ColorbarVisible','off');
% title('Potential synapses (10um)');

load ('../matFiles/ConnMatrix_CO_top500_2blocks_gamma038_08062020.mat');

trueDiag = sum(sum(ConnMat_CO_top500_2blocks_gamma038_08062020(1:size(block1_cleaned.cellIDs,2),1:size(block1_cleaned.cellIDs,2))))+ ...
    sum(sum(ConnMat_CO_top500_2blocks_gamma038_08062020(size(block1_cleaned.cellIDs,2)+1:540,size(block1_cleaned.cellIDs,2)+1:540)));

trueOffDiag = sum(sum(ConnMat_CO_top500_2blocks_gamma038_08062020(1:size(block1_cleaned.cellIDs,2),size(block1_cleaned.cellIDs,2)+1:540)))+ ...
    sum(sum(ConnMat_CO_top500_2blocks_gamma038_08062020(size(block1_cleaned.cellIDs,2)+1:540,1:size(block1_cleaned.cellIDs,2))));

%trueDiag = sum(block1_cleaned.actualSynapses11)+sum(block2_cleaned.actualSynapses22);
trueOffDiag = sum(block2_cleaned.actualSynapses21)+sum(block1_cleaned.actualSynapses12);

trueDiag_2 = sum(block1_cleaned.potentialSynpses11_2)+sum(block2_cleaned.potentialSynpses22_2);
trueOffDiag_2 = sum(block2_cleaned.potentialSynpses21_2)+sum(block1_cleaned.potentialSynpses12_2);

trueDiag_5 = sum(block1_cleaned.potentialSynpses11_5)+sum(block2_cleaned.potentialSynpses22_5);
trueOffDiag_5 = sum(block2_cleaned.potentialSynpses21_5)+sum(block1_cleaned.potentialSynpses12_5);

trueDiag_10 = sum(block1_cleaned.potentialSynpses11_10)+sum(block2_cleaned.potentialSynpses22_10);
trueOffDiag_10 = sum(block2_cleaned.potentialSynpses21_10)+sum(block1_cleaned.potentialSynpses12_10);

figure;
subplot(4,4,1)
plot([1,2,3,4],[trueDiag/trueOffDiag,trueDiag_2/trueOffDiag_2,trueDiag_5/trueOffDiag_5,trueDiag_10/trueOffDiag_10],'-ko','LineWidth',2);
set(gca,'XTick',[1:4],'XTickLabels',[{'true'},2,5,10]);
box off;
daspect([1,1,1]);
offsetAxes(gca);
xlabel('radius (\mum)');
ylabel('diagonal sum/off-diagonal sum')

% subplot(4,4,5)
% % scatter(block1_cleaned.actualSynapses11,block1_cleaned.potentialSynpses11_5);
% hold on
% % scatter(block1_cleaned.actualSynapses12,block1_cleaned.potentialSynpses12_5);
% % scatter(block2_cleaned_gamma038.actualSynapses22,block2_cleaned_gamma038.potentialSynpses22_5);
% % scatter(block2_cleaned_gamma038.actualSynapses21,block2_cleaned_gamma038.potentialSynpses21_5);
% %line([0,100],[0,100],'color','k');
% showfit(ezfit(block1_cleaned.actualSynapses11,block1_cleaned.potentialSynpses11_5,'affine'),'fitcolor',col1,'dispeqboxmode','off');
% showfit(ezfit(block1_cleaned.actualSynapses12,block1_cleaned.potentialSynpses12_5,'affine'),'fitcolor',col1,'dispeqboxmode','off','fitlinestyle',':');
% showfit(ezfit(block2_cleaned_gamma038.actualSynapses22,block2_cleaned_gamma038.potentialSynpses22_5,'affine'),'fitcolor',col2,'dispeqboxmode','off');
% showfit(ezfit(block2_cleaned_gamma038.actualSynapses21,block2_cleaned_gamma038.potentialSynpses21_5,'affine'),'fitcolor',col2,'dispeqboxmode','off','fitlinestyle',':');
% set(gca,'XLim',[0,100],'YLim',[0,100000]);
% legend({'b1-->b1','b1-->b2','b2-->b2','b2-->b1'});
% %axis square;
% xlabel('True Synapses');
% ylabel('Potential Synapses');

% peters analysis for block1_cleaned-->ABD and block2_cleaned_gamma038 -->ABD

ABDr_CellIDs = [82140,82145,82143,77648,82146,77302,77301,82194,77705,77710,77305,77672,77300,77709,77661];
ABDc_CellIDs = [82213,77654,77646,82212,77295,81172,77657,82197,77682,77154, 77652,77658,77628,77292,77688,82195,77296];
ABDIr_CellIDs = [78553,77150,77886,78547,77631,77158,77665,78552,77668,77618,77634,78556];
ABDIc_CellIDs = [78574,79051,77148,77625,77144,77692,77641,79066,77640,77643];

% get all axonal sites for ABDm neruons

motorNeurons = [ABDr_CellIDs,ABDc_CellIDs];
motorNeurons_dendPartners = [];
for i = 1:length(motorNeurons)
    [temp,prePSDid] = SynapticPartners(motorNeurons(i),1,df_cleaned);
    prePSDid = prePSDid(temp~=prePSDid(i)); % remove self touches
    motorNeurons_dendPartners = [motorNeurons_dendPartners;PrePartnerCoordinates(prePSDid,df_cleaned)];
    clear temp;
    clear prePSDid;
end
motorNeurons_dendPartners = TransformPoints(motorNeurons_dendPartners,0);

% interneuron dendrite partner sites
interNeurons = [ABDIr_CellIDs,ABDIc_CellIDs];
interNeurons_dendPartners = [];
for i = 1:length(interNeurons)
    [temp,prePSDid] = SynapticPartners(interNeurons(i),1,df_cleaned);
    prePSDid = prePSDid(temp~=prePSDid(i)); % remove self touches
    interNeurons_dendPartners = [interNeurons_dendPartners;PrePartnerCoordinates(prePSDid,df_cleaned)];
    clear temp;
    clear prePSDid;
end
interNeurons_dendPartners = TransformPoints(interNeurons_dendPartners,0);



for i = 1:length(block1_cleaned.cellIDs)
    [temp,postID] = SynapticPartners(block1_cleaned.cellIDs(i),2,df_cleaned);
    postID = postID(temp~=block1_cleaned.cellIDs(i));
    temp2= PrePartnerCoordinates(postID,df_cleaned);
    temp2 = TransformPoints(temp2,0);
    
    distMatmotor = pdist2(temp2,motorNeurons_dendPartners);
    distMatinter = pdist2(temp2,interNeurons_dendPartners);
    
    block1_cleaned.actualSynapsesMotor(i) = length(find(distMatmotor == 0));
    block1_cleaned.potentialSynapsesMotor_2(i) = length(find(distMatmotor <=2));
    block1_cleaned.potentialSynapsesMotor_5(i) = length(find(distMatmotor <=5));
    block1_cleaned.potentialSynapsesMotor_10(i) = length(find(distMatmotor <=10));
    
    block1_cleaned.actualSynapsesInter(i) = length(find(distMatinter == 0));
    block1_cleaned.potentialSynapsesInter_2(i) = length(find(distMatinter <=2));
    block1_cleaned.potentialSynapsesInter_5(i) = length(find(distMatinter <=5));
    block1_cleaned.potentialSynapsesInter_10(i) = length(find(distMatinter <=10));
    
    %block1_cleaned.axonalSites = [block1_cleaned.axonalSites;temp2]
    
    clear temp
    clear temp2
    clear postID;
    clear distMatmotor
    clear distMatinter
end


for i = 1:length(block2_cleaned.cellIDs)
    [temp,postID] = SynapticPartners(block2_cleaned.cellIDs(i),2,df_cleaned);
    postID = postID(temp~=block2_cleaned.cellIDs(i));
    temp2= PrePartnerCoordinates(postID,df_cleaned);
    temp2 = TransformPoints(temp2,0);
    
    distMatmotor = pdist2(temp2,motorNeurons_dendPartners);
    distMatinter = pdist2(temp2,interNeurons_dendPartners);
    
    block2_cleaned.actualSynapsesMotor(i) = length(find(distMatmotor == 0));
    block2_cleaned.potentialSynapsesMotor_2(i) = length(find(distMatmotor <=2));
    block2_cleaned.potentialSynapsesMotor_5(i) = length(find(distMatmotor <=5));
    block2_cleaned.potentialSynapsesMotor_10(i) = length(find(distMatmotor <=10));
    
    block2_cleaned.actualSynapsesInter(i) = length(find(distMatinter == 0));
    block2_cleaned.potentialSynapsesInter_2(i) = length(find(distMatinter <=2));
    block2_cleaned.potentialSynapsesInter_5(i) = length(find(distMatinter <=5));
    block2_cleaned.potentialSynapsesInter_10(i) = length(find(distMatinter <=10));
    
    %block2_cleaned_gamma038.axonalSites = [block2_cleaned_gamma038.axonalSites;temp2]
    
    clear temp
    clear temp2
    clear postID;
    clear distMatmotor
    clear distMatinter
end


% Potential synapse, mod to ABD

mod1ABD = sum(block1_cleaned.actualSynapsesMotor)+sum(block1_cleaned.actualSynapsesInter);
mod2ABD = sum(block2_cleaned.actualSynapsesMotor)+sum(block2_cleaned.actualSynapsesInter);

mod1ABD_2 = sum(block1_cleaned.potentialSynapsesMotor_2)+sum(block1_cleaned.potentialSynapsesInter_2);
mod1ABD_5 = sum(block1_cleaned.potentialSynapsesMotor_5)+sum(block1_cleaned.potentialSynapsesInter_5);
mod1ABD_10 = sum(block1_cleaned.potentialSynapsesMotor_10)+sum(block1_cleaned.potentialSynapsesInter_10);

mod2ABD_2 = sum(block2_cleaned.potentialSynapsesMotor_2)+sum(block2_cleaned.potentialSynapsesInter_2);
mod2ABD_5 = sum(block2_cleaned.potentialSynapsesMotor_5)+sum(block2_cleaned.potentialSynapsesInter_5);
mod2ABD_10 = sum(block2_cleaned.potentialSynapsesMotor_10)+sum(block2_cleaned.potentialSynapsesInter_10);

subplot(4,4,6)
% semilogy([mod1ABD,mod1ABD_2,mod1ABD_5,mod1ABD_10],'-o','color',[1,0.5,0],'LineWidth',2);
% hold on
% semilogy([mod2ABD,mod2ABD_2,mod2ABD_5,mod2ABD_10],'-o','color',[0,0.5,1],'LineWidth',2);

% plot([1,2,3,4],[mod1ABD,mod1ABD_2,mod1ABD_5,mod1ABD_10]./[mod2ABD,mod2ABD_2,mod2ABD_5,mod2ABD_10],...
%     '-ok','LineWidth',2);
% yyaxis right
plot([1,2,3,4],[mod2ABD,mod2ABD_2,mod2ABD_5,mod2ABD_10]./[mod1ABD,mod1ABD_2,mod1ABD_5,mod1ABD_10],...
    '-ok','LineWidth',2);
%daspect([1,2,1]);
box off;
set(gca,'XTick',[1:4],'XTickLabels',[{'true'},2,5,10],'YLIM',[0,20]);
ylabel('mod1/mod2 synapses');
xlabel('radius (\mum)');
daspect([1,3,1]);
offsetAxes(gca);

%axis square;
%legend({'mod1','mod2'},'Location','bestoutside');

%% calculate recurrent fractions

for i = 1:length(block2_cleaned.cellIDs)
    A = SynapticPartners(block2_cleaned.cellIDs(i),1,df_cleaned);
    block2_cleaned.recurrentFraction(i) = sum(ismember(A,block2_cleaned.cellIDs))./length(A);
    clear A;
end

for i = 1:length(block1_cleaned.cellIDs)
    A = SynapticPartners(block1_cleaned.cellIDs(i),1,df_cleaned);
    block1_cleaned.recurrentFraction(i) = sum(ismember(A,block1_cleaned.cellIDs))./length(A);
    clear A;
end

% submodules
for i = 1:length(block1_cleaned_submod.cellIDs)
    A = SynapticPartners(block1_cleaned_submod.cellIDs(i),1,df_cleaned);
    block1_cleaned_submod.recurrentFraction(i) = sum(ismember(A,block1_cleaned_submod.cellIDs))./length(A);
    clear A;
end

for i = 1:length(block2_cleaned_submod.cellIDs)
    A = SynapticPartners(block2_cleaned_submod.cellIDs(i),1,df_cleaned);
    block2_cleaned_submod.recurrentFraction(i) = sum(ismember(A,block2_cleaned_submod.cellIDs))./length(A);
    clear A;
end


% feedf_cleanedorward 
%%
figure;
subplot(4,4,5)
%histogram(rmoutliers(block2_cleaned.recurrentFraction),'BinWidth',0.025,'EdgeColor',[0,0.5,1],'DisplayStyle','stairs','LineWidth',2);
histogram(block2_cleaned.recurrentFraction,'BinWidth',0.025,'EdgeColor',[0,0.5,1],'DisplayStyle','stairs','LineWidth',2);
hold on;
%histogram(rmoutliers(block1_cleaned.recurrentFraction),'BinWidth',0.025,'EdgeColor',[1,0.5,0],'DisplayStyle','stairs','LineWidth',2);
histogram(block1_cleaned.recurrentFraction,'BinWidth',0.025,'EdgeColor',[1,0.5,0],'DisplayStyle','stairs','LineWidth',2);
box off;
legend({'mod2','mod1'},'Location','bestoutside');
%axis square;
xlabel(' Recurrent fraction');
ylabel(' count');
%offsetAxes(gca);
[h,p] = kstest2(block2_cleaned.recurrentFraction,block1_cleaned.recurrentFraction);
%title(p)

subplot(4,4,1)
grp1 = ones(size(block1_cleaned.recurrentFraction,2),1);
grp2 = 2*ones(size(block2_cleaned.recurrentFraction,2),1);
boxplot([block1_cleaned.recurrentFraction,block2_cleaned.recurrentFraction],[grp1;grp2],...
    'Orientation','horizontal','Colors',[1,0.5,0;0,0.5,1]);
box off;
daspect([1,20,1])


block1_cleaned_submod.color = [0.6207    0.3103    0.2759];
block2_cleaned_submod.color = [0.5172    0.5172    1.0000];

subplot(4,4,3)
grp1 = ones(size(block1_cleaned_submod.recurrentFraction,2),1);
grp2 = 2*ones(size(block2_cleaned_submod.recurrentFraction,2),1);
boxplot([block1_cleaned_submod.recurrentFraction,block2_cleaned_submod.recurrentFraction],[grp1;grp2],...
    'Orientation','horizontal','Colors',[block1_cleaned_submod.color;block2_cleaned_submod.color]);
box off;
daspect([1,20,1])


subplot(4,4,7)
histogram(block1_cleaned_submod.recurrentFraction,'BinWidth',0.025,'EdgeColor',block1_cleaned_submod.color,'DisplayStyle','stairs','LineWidth',2);
hold on;
histogram(block2_cleaned_submod.recurrentFraction,'BinWidth',0.025,'EdgeColor',block2_cleaned_submod.color,'DisplayStyle','stairs','LineWidth',2);
legend({'mod2a','mod2b'},'Location','bestoutside');
xlabel(' Recurrent fraction');
ylabel(' count');
box off;
%offsetAxes(gca);
%[h,p] = kstest2(block1_cleaned_submod.recurrentFraction,block2_cleaned_submod.recurrentFraction);
%title(p)

%%
%cbl = [80478,77251,80529,80341,78236,80511,80482,80502,80552,80514,80554,80352,80506,80496,80519,80490,80518,80339,78239,80524,78280,80547];

% 
% matOrder_commdet = [DOs,...
%     block1_cleaned,block2_cleaned_gamma038,Block3,...
%     ABDr_CellIDs,ABDc_CellIDs,ABDIr_CellIDs,ABDIc_CellIDs];
%    
% [~,matIndex_commdet] = ismember(matOrder_commdet,AllCells);
% 
% connMatForPlot = ConnMatrixPre(matIndex_commdet,matIndex_commdet);
% 
% for i = 1:size(matOrder_commdet,2)
%     matOrder_commdet_totalInputs(i) = length(SynapticPartners(matOrder_commdet(i),1,df_cleaned));
%     connMatForPlot_Norm(i,:) = connMatForPlot(i,:)./matOrder_commdet_totalInputs(i);
% end
% connMatForPlot_Norm(connMatForPlot_Norm>0.05) = 0.05;

%% plot loation of neurons in block5 and block6

% transform_swc_AV(block5_gamma09.cellIDss,block5_gamma09.color,[],true,false);
% transform_swc_AV(block6_gamma09.cellIDss,block6_gamma09.color,[],true,false);
% figure;
% transform_swc_AV(block1_cleaned_submod.cellIDs,hex2rgb('#cc5643'),[],true,false);
% transform_swc_AV(block2_cleaned_submod.cellIDs,hex2rgb('#8773c9'),[],false,false);

block2_cleaned.removeSoma = [78617,81066,78560,81329,80990,81021];
block2_cleaned.noSoma = [78583,78883,77122,78406,79746,77581,76690,80548,77848,78860,77447,77872,79072,80606,77651,81374,81817,78443,78544,79055,77163,81297,77132,77433,81407,77349,80974,78351,77140,79062,77344,81537,78441,76697,81391,78629,77467,78150,78650,78226,80629,80995,76677,81431,76673,77797,79042,80539,76618,78421,79074,78297,77607,77162,79069,79022,78858,81793,80647,81317,81060,81007,79046,80757,77805,80743,80472,76701,80294,78641,79743,78649,81765,78254,77121,81322,77151,79720,79976,77329,77868,78540,80510,76622,76661,78601,81792,78543,81161,76627,81293,78246,77146,78911,78255,81774,78356,80800,81423,77434,77465,80679,77435,77341,79054,78357,77592,78949,77845,81336,81311,77591,77239,78667,78157,77124,80956,79722,77327,77844,77621,81637,81400,78646,79085,81002,77460,79080,77437,79048,77142,77389,78241,80947,78404,77636,77821,81312,79953,80681,79058,79033,77822,80746,80750,77656,80850,79067,80625,78453,76623,81363,81550,80626,81683,78633,79044,81395,77152,81443,80262,80315,80821,77602,77251,81408,77369,81580,81417,76625,78696,77580,79852,77126,81295,81559,81338];
block2_cleaned.cleanSoma = block2_cleaned.cellIDs(~ismember(block2_cleaned.cellIDs,block2_cleaned.removeSoma));

% get origins
block1_cleaned_submod.cleanedSoma = block1_cleaned_submod.cellIDs(ismember(block1_cleaned_submod.cellIDs,block2_cleaned.cleanSoma));
block2_cleaned_submod.cleanedSoma = block2_cleaned_submod.cellIDs(ismember(block2_cleaned_submod.cellIDs,block2_cleaned.cleanSoma));

block1_cleaned_submod.noSomaIndex = ismember(block1_cleaned_submod.cleanedSoma,block2_cleaned.noSoma);
block2_cleaned_submod.noSomaIndex = ismember(block2_cleaned_submod.cleanedSoma,block2_cleaned.noSoma);

block1_cleaned_submod.inVolIndex = ~ismember(block1_cleaned_submod.cleanedSoma,block2_cleaned.noSoma);
block2_cleaned_submod.inVolIndex = ~ismember(block2_cleaned_submod.cleanedSoma,block2_cleaned.noSoma);

block1_cleaned_submod.cleanOrigins = getOrigin(block1_cleaned_submod.cleanedSoma);
block2_cleaned_submod.cleanOrigins = getOrigin(block2_cleaned_submod.cleanedSoma);

% scatter histogram of spatial locations


figure;
h = scatterhist([block1_cleaned_submod.cleanOrigins(block1_cleaned_submod.inVolIndex,1);block2_cleaned_submod.cleanOrigins(block2_cleaned_submod.inVolIndex,1)],...
    [block1_cleaned_submod.cleanOrigins(block1_cleaned_submod.inVolIndex,2);block2_cleaned_submod.cleanOrigins(block2_cleaned_submod.inVolIndex,2)],...
    'Group',[ones(sum(block2_cleaned_submod.inVolIndex),1);2*ones(sum(block1_cleaned_submod.inVolIndex),1)],...
'Kernel','on','color',[block1_cleaned_submod.color;block2_cleaned_submod.color],...
    'Marker','o','MarkerSize',5,'LineWidth',2);
set(h(1),'YDir','reverse','YLim',[400,800],'XLim',[250,400]);
set(h(3),'XDir','reverse','Visible','on');
set(h(2),'Visible','on');
daspect([1,1,1]);
hold on;
scatter(h(1),block1_cleaned_submod.cleanOrigins(block1_cleaned_submod.noSomaIndex,1),...
    block1_cleaned_submod.cleanOrigins(block1_cleaned_submod.noSomaIndex,2),...
    20,'MarkerEdgeColor',block1_cleaned_submod.color,'Marker','o','LineWidth',0.25);
scatter(h(1),block2_cleaned_submod.cleanOrigins(block2_cleaned_submod.noSomaIndex,1),...
    block2_cleaned_submod.cleanOrigins(block2_cleaned_submod.noSomaIndex,2),...
    20,'MarkerEdgeColor',block2_cleaned_submod.color,'Marker','o','LineWidth',0.25);
daspect([1,1,1]);



% location of presynapsess

block1_cleaned_submod.dendSites = [];
block1_cleaned_submod.dendpartnerSites = [];
for i = 1:length(block1_cleaned_submod.cellIDs)
    [temp,preIDs] = SynapticPartners(block1_cleaned_submod.cellIDs(i),1,df_cleaned);
    preIDs = preIDs(temp~=block1_cleaned_submod.cellIDs(i)); % remove self touches
    block1_cleaned_submod.dendSites = [block1_cleaned_submod.dendSites;PrePartnerCoordinates(preIDs,df_cleaned)];
    block1_cleaned_submod.dendpartnerSites = [block1_cleaned_submod.dendpartnerSites;PostPartnerCoordinates(preIDs,df_cleaned)];
    clear preIDs;
    clear temp;
end
block1_cleaned_submod.dendSites = TransformPoints(block1_cleaned_submod.dendSites,0);
block1_cleaned_submod.dendpartnerSites = TransformPoints(block1_cleaned_submod.dendpartnerSites,0);


block2_cleaned_submod.dendSites = [];
block2_cleaned_submod.dendpartnerSites = [];
for i = 1:length(block2_cleaned_submod.cellIDs)
    [temp,preIDs] = SynapticPartners(block2_cleaned_submod.cellIDs(i),1,df_cleaned);
    preIDs = preIDs(temp~=block2_cleaned_submod.cellIDs(i)); % remove self touches
    block2_cleaned_submod.dendSites = [block2_cleaned_submod.dendSites;PrePartnerCoordinates(preIDs,df_cleaned)];
    block2_cleaned_submod.dendpartnerSites = [block2_cleaned_submod.dendpartnerSites;PostPartnerCoordinates(preIDs,df_cleaned)];
    clear preIDs;
    clear temp;
end
block2_cleaned_submod.dendSites = TransformPoints(block2_cleaned_submod.dendSites,0);
block2_cleaned_submod.dendpartnerSites = TransformPoints(block2_cleaned_submod.dendpartnerSites,0);


% location of postsynapses
block1_cleaned_submod.axonSites = [];
block1_cleaned_submod.axonPartnerSites = [];
for i = 1:length(block1_cleaned_submod.cellIDs)
    [temp,preIDs] = SynapticPartners(block1_cleaned_submod.cellIDs(i),2,df_cleaned);
    preIDs = preIDs(temp~=block1_cleaned_submod.cellIDs(i)); % remove self touches
    block1_cleaned_submod.axonSites = [block1_cleaned_submod.axonSites;PostPartnerCoordinates(preIDs,df_cleaned)];
    block1_cleaned_submod.axonPartnerSites = [block1_cleaned_submod.axonPartnerSites;PrePartnerCoordinates(preIDs,df_cleaned)];
    clear preIDs;
end
block1_cleaned_submod.axonSites = TransformPoints(block1_cleaned_submod.axonSites,0);
block1_cleaned_submod.axonPartnerSites = TransformPoints(block1_cleaned_submod.axonPartnerSites,0);

block2_cleaned_submod.axonSites = [];
block2_cleaned_submod.axonPartnerSites = [];
for i = 1:length(block2_cleaned_submod.cellIDs)
    [temp,preIDs] = SynapticPartners(block2_cleaned_submod.cellIDs(i),2,df_cleaned);
    preIDs = preIDs(temp~=block2_cleaned_submod.cellIDs(i)); % remove self touches
    block2_cleaned_submod.axonSites = [block2_cleaned_submod.axonSites;PostPartnerCoordinates(preIDs,df_cleaned)];
    block2_cleaned_submod.axonPartnerSites = [block2_cleaned_submod.axonPartnerSites;PrePartnerCoordinates(preIDs,df_cleaned)];
    clear preIDs;
end
block2_cleaned_submod.axonSites = TransformPoints(block2_cleaned_submod.axonSites,0);
block2_cleaned_submod.axonPartnerSites = TransformPoints(block2_cleaned_submod.axonPartnerSites,0);


% plot histograms of the inputs
Xcords_Pre = [block1_cleaned_submod.dendSites(1:5:end,1);block2_cleaned_submod.dendSites(1:5:end,1)];
Ycords_Pre = [block1_cleaned_submod.dendSites(1:5:end,2);block2_cleaned_submod.dendSites(1:5:end,2)];
Zcords_pre = [block1_cleaned_submod.dendSites(1:5:end,3);block2_cleaned_submod.dendSites(1:5:end,3)];

GpIDs_Pre = [ones(length(block1_cleaned_submod.dendSites(1:5:end,1)),1);2*ones(length(block2_cleaned_submod.dendSites(1:5:end,1)),1)];

figure;
h = scatterhist(Xcords_Pre,Ycords_Pre,'Group',GpIDs_Pre,'Kernel','on','color',[block1_cleaned_submod.color;block2_cleaned_submod.color],...
    'Marker','.','MarkerSize',1,'LineWidth',2);
set(h(1),'YDir','reverse','YLim',[400,800],'XLim',[250,400]);
set(h(3),'XDir','reverse');
daspect([1,1,1]);

% plot histogram of the outputs
Xcords_Post = [block1_cleaned_submod.axonSites(1:5:end,1);block2_cleaned_submod.axonSites(1:5:end,1)];
Ycords_Post = [block1_cleaned_submod.axonSites(1:5:end,2);block2_cleaned_submod.axonSites(1:5:end,2)];
Zcords_Post = [block1_cleaned_submod.axonSites(1:5:end,3);block2_cleaned_submod.axonSites(1:5:end,3)];
GpIDs_Post = [ones(length(block1_cleaned_submod.axonSites(1:5:end,1)),1);2*ones(length(block2_cleaned_submod.axonSites(1:5:end,1)),1)];

figure;
h = scatterhist(Xcords_Post,Ycords_Post,'Group',GpIDs_Post,'Kernel','on','color',[block1_cleaned_submod.color;block2_cleaned_submod.color],...
    'Marker','.','MarkerSize',0.5,'LineWidth',2);
set(h(1),'YDir','reverse','YLim',[400,800],'XLim',[250,400]);
set(h(3),'XDir','reverse');
daspect([1,1,1]);

%% Peters analysis for submodules; block5 -->ABDm, block6 -->ABDi

% block1_ABDi --> block1_ABDi
for i = 1:length(block1_cleaned_submod.cellIDs)
    [prePartner,prePartnerPSD] = SynapticPartners(block1_cleaned_submod.cellIDs(i),1,df_cleaned);
    % remove self touches
    %block1_cleaned.block1_cleanedtoblock1_cleanedActual(i) = sum(ismember(prePartner,setdiff(block1_cleaned.cellIDss,block1_cleaned.cellIDss(i)))); 
    prePartnerPSD = prePartnerPSD(prePartner~=block1_cleaned_submod.cellIDs(i));

    inputs = PrePartnerCoordinates(prePartnerPSD,df_cleaned); % location of presynaptic partner, on the axon
    inputs = TransformPoints(inputs,0);
    % distance matrix block1_cleaned -->block1_cleaned
    distMat11 = pdist2(block1_cleaned_submod.axonPartnerSites,inputs);
    block1_cleaned_submod.actualSynapses11(i) = length(find(distMat11 == 0));
    block1_cleaned_submod.potentialSynpses11_05(i) = length(find(distMat11 <= 0.5));
    block1_cleaned_submod.potentialSynpses11_1(i) = length(find(distMat11 <= 1));
    block1_cleaned_submod.potentialSynpses11_2(i) = length(find(distMat11 <= 2));
    block1_cleaned_submod.potentialSynpses11_5(i) = length(find(distMat11 <= 5));
    block1_cleaned_submod.potentialSynpses11_10(i) = length(find(distMat11 <=10));


    %distance Matrix block1_cleaned -->block2_cleaned_gamma038
    distMat12 = pdist2(block2_cleaned_submod.axonPartnerSites,inputs);
    block1_cleaned_submod.actualSynapses12(i) = length(find(distMat12 == 0));
    block1_cleaned_submod.potentialSynpses12_05(i) = length(find(distMat12 <= 0.5));
    block1_cleaned_submod.potentialSynpses12_1(i) = length(find(distMat12 <= 1));
    block1_cleaned_submod.potentialSynpses12_2(i) = length(find(distMat12 <= 2));
    block1_cleaned_submod.potentialSynpses12_5(i) = length(find(distMat12 <= 5));
    block1_cleaned_submod.potentialSynpses12_10(i) = length(find(distMat12 <= 10));

    
    clear PrePartner;
    clear inputs;
    clear distMat11;
    clear distMat12;
end

% block2_ABDm --> block2_ABDm
for i = 1:length(block2_cleaned_submod.cellIDs)
    [prePartner,prePartnerPSD] = SynapticPartners(block2_cleaned_submod.cellIDs(i),1,df_cleaned);
    % remove self touches
    %block1_cleaned.block1_cleanedtoblock1_cleanedActual(i) = sum(ismember(prePartner,setdiff(block1_cleaned.cellIDss,block1_cleaned.cellIDss(i)))); 
    prePartnerPSD = prePartnerPSD(prePartner~=block2_cleaned_submod.cellIDs(i));

    inputs = PrePartnerCoordinates(prePartnerPSD,df_cleaned); % location of presynaptic partner, on the axon
    inputs = TransformPoints(inputs,0);
    % distance matrix block1_cleaned -->block1_cleaned
    distMat22 = pdist2(block2_cleaned_submod.axonPartnerSites,inputs);
    block2_cleaned_submod.actualSynapses22(i) = length(find(distMat22 == 0));
    block2_cleaned_submod.potentialSynpses22_05(i) = length(find(distMat22 <= 0.5));
    block2_cleaned_submod.potentialSynpses22_1(i) = length(find(distMat22 <= 1));
    block2_cleaned_submod.potentialSynpses22_2(i) = length(find(distMat22 <= 2));
    block2_cleaned_submod.potentialSynpses22_5(i) = length(find(distMat22 <= 5));
    block2_cleaned_submod.potentialSynpses22_10(i) = length(find(distMat22 <=10));


    %distance Matrix block1_cleaned -->block2_cleaned_gamma038
    distMat21 = pdist2(block1_cleaned_submod.axonPartnerSites,inputs);
    block2_cleaned_submod.actualSynapses21(i) = length(find(distMat21 == 0));
    block2_cleaned_submod.potentialSynpses21_05(i) = length(find(distMat21 <= 0.5));
    block2_cleaned_submod.potentialSynpses21_1(i) = length(find(distMat21 <= 1));
    block2_cleaned_submod.potentialSynpses21_2(i) = length(find(distMat21 <= 2));
    block2_cleaned_submod.potentialSynpses21_5(i) = length(find(distMat21 <= 5));
    block2_cleaned_submod.potentialSynpses21_10(i) = length(find(distMat21 <= 10));

    
    clear PrePartner;
    clear inputs;
    clear distMat22;
    clear distMat21;
end

%block1_cleaned_submod.indices = 

load ConnMatrix_CO_OM_top500_2blocks_gamma038_08072020.mat

trueDiag_submods = sum(sum(ConnMat_CO_OM_top500_2blocks_gamma038_08072020(1:size(block1_cleaned_submod.cellIDs,2),1:size(block1_cleaned_submod.cellIDs,2))))+ ...
    sum(sum(ConnMat_CO_OM_top500_2blocks_gamma038_08072020(size(block1_cleaned_submod.cellIDs,2)+1:289,size(block1_cleaned_submod.cellIDs,2)+1:289)));

trueOffDiag_submods = sum(sum(ConnMat_CO_OM_top500_2blocks_gamma038_08072020(1:size(block1_cleaned_submod.cellIDs,2),size(block1_cleaned_submod.cellIDs,2)+1:289)))+ ...
    sum(sum(ConnMat_CO_OM_top500_2blocks_gamma038_08072020(size(block1_cleaned_submod.cellIDs,2)+1:289,1:size(block1_cleaned_submod.cellIDs,2))));


% trueDiag_submods = sum(block1_cleaned_submod.actualSynapses11)+sum(block2_cleaned_submod.actualSynapses22);
% trueOffDiag_submods = sum(block1_cleaned_submod.actualSynapses12) + sum( block2_cleaned_submod.actualSynapses21);

trueDiag_submods_2 = sum(block1_cleaned_submod.potentialSynpses11_2)+sum(block2_cleaned_submod.potentialSynpses22_2);
trueOffDiag_submods_2 = sum(block1_cleaned_submod.potentialSynpses12_2) + sum( block2_cleaned_submod.potentialSynpses21_2);

trueDiag_submods_5 = sum(block1_cleaned_submod.potentialSynpses11_5)+sum(block2_cleaned_submod.potentialSynpses22_5);
trueOffDiag_submods_5 = sum(block1_cleaned_submod.potentialSynpses12_5) + sum( block2_cleaned_submod.potentialSynpses21_5);

trueDiag_submods_10 = sum(block1_cleaned_submod.potentialSynpses11_10)+sum(block2_cleaned_submod.potentialSynpses22_10);
trueOffDiag_submods_10 = sum(block1_cleaned_submod.potentialSynpses12_10) + sum( block2_cleaned_submod.potentialSynpses21_10);


figure;

subplot(4,4,1)

plot([1,2,3,4],[trueDiag_submods/trueOffDiag_submods,trueDiag_submods_2/trueOffDiag_submods_2,...
    trueDiag_submods_5/trueOffDiag_submods_5,trueDiag_submods_10/trueOffDiag_submods_10],'-ok','LineWidth',2);
set(gca,'XTick',[1:4],'XTickLabels',[{'true'},2,5,10],'YLim',[1,6]);
box off;
daspect([1,1,1]);
xlabel('radius (\mum)');
ylabel('$\displaystyle\frac{sum(Diag.)}{sum(off-Diag.)}$','Interpreter','latex')

% get all block1_ABDi-->motor

block1_cleaned_submod.axonalSites = [];
for i = 1:length(block1_cleaned_submod.cellIDs)
    [temp,postID] = SynapticPartners(block1_cleaned_submod.cellIDs(i),2,df_cleaned);
    postID = postID(temp~=block1_cleaned_submod.cellIDs(i));
    temp2= PrePartnerCoordinates(postID,df_cleaned);
    temp2 = TransformPoints(temp2,0);
    
    distMatmotor = pdist2(temp2,motorNeurons_dendPartners);
    distMatinter = pdist2(temp2,interNeurons_dendPartners);
    
    block1_cleaned_submod.actualSynapsesMotor(i) = length(find(distMatmotor == 0));
    block1_cleaned_submod.potentialSynapsesMotor_2(i) = length(find(distMatmotor <=2));
    block1_cleaned_submod.potentialSynapsesMotor_5(i) = length(find(distMatmotor <=5));
    block1_cleaned_submod.potentialSynapsesMotor_10(i) = length(find(distMatmotor <=10));
    
    block1_cleaned_submod.actualSynapsesInter(i) = length(find(distMatinter == 0));
    block1_cleaned_submod.potentialSynapsesInter_2(i) = length(find(distMatinter <=2));
    block1_cleaned_submod.potentialSynapsesInter_5(i) = length(find(distMatinter <=5));
    block1_cleaned_submod.potentialSynapsesInter_10(i) = length(find(distMatinter <=10));
    
    
    block1_cleaned_submod.axonalSites = [block1_cleaned_submod.axonalSites;temp2];
    
    clear temp
    clear postID;
end


% block6 --> motor

block2_cleaned_submod.axonalSites = [];
for i = 1:length(block2_cleaned_submod.cellIDs)
    [temp,postID] = SynapticPartners(block2_cleaned_submod.cellIDs(i),2,df_cleaned);
    postID = postID(temp~=block2_cleaned_submod.cellIDs(i));
    temp2= PrePartnerCoordinates(postID,df_cleaned);
    temp2 = TransformPoints(temp2,0);
    
    distMatmotor = pdist2(temp2,motorNeurons_dendPartners);
    distMatinter = pdist2(temp2,interNeurons_dendPartners);
    
    block2_cleaned_submod.actualSynapsesMotor(i) = length(find(distMatmotor == 0));
    block2_cleaned_submod.potentialSynapsesMotor_2(i) = length(find(distMatmotor <=2));
    block2_cleaned_submod.potentialSynapsesMotor_5(i) = length(find(distMatmotor <=5));
    block2_cleaned_submod.potentialSynapsesMotor_10(i) = length(find(distMatmotor <=10));
    
    block2_cleaned_submod.actualSynapsesInter(i) = length(find(distMatinter == 0));
    block2_cleaned_submod.potentialSynapsesInter_2(i) = length(find(distMatinter <=2));
    block2_cleaned_submod.potentialSynapsesInter_5(i) = length(find(distMatinter <=5));
    block2_cleaned_submod.potentialSynapsesInter_10(i) = length(find(distMatinter <=10));
    
    
    block2_cleaned_submod.axonalSites = [block2_cleaned_submod.axonalSites;temp2];
    
    
    clear temp
    clear postID;
end


subplot(4,4,2)
plot([1,2,3,4],[sum(block2_cleaned_submod.actualSynapsesMotor),sum(block2_cleaned_submod.potentialSynapsesMotor_2),...
    sum(block2_cleaned_submod.potentialSynapsesMotor_5),sum(block2_cleaned_submod.potentialSynapsesMotor_10)]./...
    [sum(block1_cleaned_submod.actualSynapsesMotor),sum(block1_cleaned_submod.potentialSynapsesMotor_2),...
    sum(block1_cleaned_submod.potentialSynapsesMotor_5),sum(block1_cleaned_submod.potentialSynapsesMotor_10)],...
    '-ok','LineWidth',2);
set(gca,'XTick',[1:4],'XTickLabels',[{'true'},2,5,10]);
%ylabel('mod2b/mod2a on ABD_M');
box off;
offsetAxes(gca);
hold on;
%subplot(4,4,3)
plot([1,2,3,4],[sum(block1_cleaned_submod.actualSynapsesInter),sum(block1_cleaned_submod.potentialSynapsesInter_2),...
    sum(block1_cleaned_submod.potentialSynapsesInter_5),sum(block1_cleaned_submod.potentialSynapsesInter_10)]./...
    [sum(block2_cleaned_submod.actualSynapsesInter),sum(block2_cleaned_submod.potentialSynapsesInter_2),...
    sum(block2_cleaned_submod.potentialSynapsesInter_5),sum(block2_cleaned_submod.potentialSynapsesInter_10)],...
    '-or','LineWidth',2);
%set(gca,'XTick',[1:4],'XTickLabels',[{'true'},2,5,10]);
%ylabel('mod2a/mod2b onto ABD_I');
legend({'mod2b/mod2a -> ABD_M','mod2a/mod2b -> ABD_I'})
box off;
%offsetAxes(gca);

subplot(4,4,3)
plot([1,2,3,4],[(sum(block1_cleaned_submod.actualSynapsesInter)+sum(block2_cleaned_submod.actualSynapsesMotor))/(sum(block1_cleaned_submod.actualSynapsesMotor)+sum( block2_cleaned_submod.actualSynapsesInter));...
                (sum(block2_cleaned_submod.potentialSynapsesInter_2)+sum(block2_cleaned_submod.potentialSynapsesMotor_2))/(sum(block1_cleaned_submod.potentialSynapsesMotor_2)+sum(block2_cleaned_submod.potentialSynapsesInter_2));...
                (sum(block2_cleaned_submod.potentialSynapsesInter_5)+sum(block2_cleaned_submod.potentialSynapsesMotor_5))/(sum(block1_cleaned_submod.potentialSynapsesMotor_5)+sum(block2_cleaned_submod.potentialSynapsesInter_5));...
                (sum(block2_cleaned_submod.potentialSynapsesInter_10)+sum(block2_cleaned_submod.potentialSynapsesMotor_10))/(sum(block1_cleaned_submod.potentialSynapsesMotor_10)+sum(block2_cleaned_submod.potentialSynapsesInter_10))],...
                '-ko','LineWidth',2);
set(gca,'XTick',[1:4],'XTickLabels',[{'true'},2,5,10],'YLim',[1,10]);
ylabel('sum(off-Diag.)/sum(Diag)');
box off;
%daspect([1,2,1])
%offsetAxes(gca);

subplot(4,4,4)
plot([1,2,3,4],[(sum(block1_cleaned_submod.actualSynapsesMotor)+sum( block2_cleaned_submod.actualSynapsesInter))/(sum(block1_cleaned_submod.actualSynapsesInter)+sum(block2_cleaned_submod.actualSynapsesMotor));...
                (sum(block1_cleaned_submod.potentialSynapsesMotor_2)+sum(block2_cleaned_submod.potentialSynapsesInter_2))/(sum(block2_cleaned_submod.potentialSynapsesInter_2)+sum(block2_cleaned_submod.potentialSynapsesMotor_2));...
                (sum(block1_cleaned_submod.potentialSynapsesMotor_5)+sum(block2_cleaned_submod.potentialSynapsesInter_5))/(sum(block2_cleaned_submod.potentialSynapsesInter_5)+sum(block2_cleaned_submod.potentialSynapsesMotor_5));...
                (sum(block1_cleaned_submod.potentialSynapsesMotor_10)+sum(block2_cleaned_submod.potentialSynapsesInter_10))/(sum(block2_cleaned_submod.potentialSynapsesInter_10)+sum(block2_cleaned_submod.potentialSynapsesMotor_10))],...
                '-ko','LineWidth',2);
set(gca,'XTick',[1:4],'XTickLabels',[{'true'},2,5,10]);
ylabel('sum(Diag.)/sum(off-Diag)');
box off;
%daspect([10,1,1])

% figure;
% subplot(4,4,1)
% heatmap([mean(block1_cleaned_submod.actualSynapsesMotor),mean(block1_cleaned_submod.actualSynapsesInter);mean(block2_cleaned_submod.actualSynapsesMotor),mean( block2_cleaned_submod.actualSynapsesInter)],...
%     'ColorbarVisible','off','XData',{'ABDm','ABDi'},'YData',{'Int_ABDi','Int_ABDm'});
% title('True synapses');
% 
% subplot(4,4,2)
% heatmap([mean(block1_cleaned_submod.potentialSynapsesMotor_5),mean(block1_cleaned_submod.potentialSynapsesInter_5);mean(block2_cleaned_submod.potentialSynapsesMotor_5),mean(block2_cleaned_submod.potentialSynapsesInter_5)],...
%     'ColorbarVisible','off','XData',{'ABDm','ABDi'},'YData',{'Int_ABDi','Int_ABDm'});
% title('Potential synapses (5um)');


%% plot RS cell locations in block3, block4

% Rov3 = [79961,81410,77449,83151,77456,77441,81611];
% MiV2 = [78577,77694,77695,78566];
% MiV1 = [80327,82221,82220,78940,82218,82217];
% MiM1 = [77265];
% RoM2 = [79244];
% Mid3i = [76562];
% RoL1 = [79395];
% 
% temp = distinguishable_colors(17,'w');
% RSall = [Rov3,MiV2,MiV1,MiM1,RoM2,Mid3i,RoL1];
% McellSoma = [66732, 25866, 16807];
% McellSoma = TransformPoints(McellSoma,0);
% transform_swc_AV(77099,temp(17,:),[],false,false);
% transform_swc_AV(block3_gamma09_WithNull.cellIDss(ismember(block3_gamma09_WithNull.cellIDss,RSall)),temp(11,:),[],false,false);
% transform_swc_AV(block4_gamma09_WithNull.cellIDss(ismember(block4_gamma09_WithNull.cellIDss,RSall)),temp(12,:),[],true,false);
% scatter3(McellSoma(1),McellSoma(2),McellSoma(3),180,'MarkerEdgeColor',temp(17,:),'MarkerFaceColor',temp(17,:));
% 
% 
% % transform_swc_AV(block3_gamma09_WithNull.cellIDss(ismember(block3_gamma09_WithNull.cellIDss,RSall)),hex2rgb('#91b4ba'),[],false,false);
% % transform_swc_AV(block4_gamma09_WithNull.cellIDss(ismember(block4_gamma09_WithNull.cellIDss,RSall)),hex2rgb('#bb9a4e'),[],true,false);
% 
% %%
% 
% for i = 1:length(block2_cleaned_gamma038.cellIDss)
%     A = SynapticPartners(block2_cleaned_gamma038.cellIDss(i),1,df_cleaned);
%     block2_cleaned_gamma038.recurrentFraction(i) = sum(ismember(A,block2_cleaned_gamma038.cellIDss))./length(A);
%     block2_cleaned_gamma038.recurrentSynapses(i) = sum(ismember(A,block2_cleaned_gamma038.cellIDss));
%     clear A;
% end
% 
% for i = 1:length(block1_cleaned.cellIDss)
%     A = SynapticPartners(block1_cleaned.cellIDss(i),1,df_cleaned);
%     block1_cleaned.recurrentFraction(i) = sum(ismember(A,block1_cleaned.cellIDss))./length(A);
%     block1_cleaned.recurrentSynapses(i)  = sum(ismember(A,block1_cleaned.cellIDss));
%     clear A;
% end
% 
% % feedf_cleanedorward
% 
% figure;
% subplot(4,4,1)
% histogram(block2_cleaned_gamma038.recurrentFraction,'FaceColor',[0,0.5,1]);
% hold on;
% histogram(block1_cleaned.recurrentFraction,'FaceColor',[1,0.5,0]);
% legend({'Int','non-Int'});
% axis square;
% xlabel('% recurrent fraction');
% box off;
% offsetAxes(gca);
% %%
% 
% [a,b,c] = OSI(block1_cleaned.cellIDss,df_cleaned);
% [d,e,f] = OSI(block2_cleaned_gamma038.cellIDss,df_cleaned);
% 
% g  = b+c;
% h  = e+f;
% 
% 
% 
% block1_cleaned.pdf_cleaned =  fitdist(g(g>0),'Kernel');
% block2_cleaned_gamma038.pdf_cleaned = fitdist(h(h>0),'Kernel');
% 
% subplot(4,4,1)
% % plot(0:1:100,pdf_cleaned(block1_cleaned.pdf_cleaned,0:1:100),'-','color',[1,0.5,0],'LineWidth',2);
% hold on;
% % yyaxis right;
% histogram(g(g>0),'BinWidth',5,'EdgeColor',[1,0.5,0],'DisplayStyle','stairs','LineWidth',2);
% % yyaxis left
% % plot(0:1:200,pdf_cleaned(block2_cleaned_gamma038.pdf_cleaned,0:1:200),'-','color',[0,0.5,1],'LineWidth',2);
% % yyaxis right
% histogram(h(h>0),'BinWidth',5,'EdgeColor',[0,0.5,1],'DisplayStyle','stairs','LineWidth',2);
% 
% axis square;
% box off;
% xlabel('ABD_M + ABD_I synapses');
% ylabel('Count');
% offsetAxes(gca);
% 
% [h1,p1] = kstest2(g(g>0),h(h>0));

%% Synapses size


% mean synapse size in block1_cleaned (i-->j)

for i= 1:length(block1_cleaned.cellIDs)
    if isExistReRoot(block1_cleaned.cellIDs(i))
        block1_cleaned.tree(i) = SwctoZbrian(block1_cleaned.cellIDs(i));
    end
end

block1_cleaned.b1tob1size = [];
block1_cleaned.b1tob1postSynCoord = [];
block1_cleaned.b1tob1postSynPathLength = [];
block1_cleaned.b1tob1postSynPathLength_norm = [];

for i = 1:length(block1_cleaned.cellIDs) % axon
    for j = 1:length(block1_cleaned.cellIDs) % dendrite
        if ~isequal(block1_cleaned.cellIDs(i),block1_cleaned.cellIDs(j))
            temp = df_cleaned.size(df_cleaned.presyn_segid==block1_cleaned.cellIDs(i) & df_cleaned.postsyn_segid == block1_cleaned.cellIDs(j));
            block1_cleaned.b1tob1size = [block1_cleaned.b1tob1size;temp];
            tempx = df_cleaned.postsyn_x(df_cleaned.presyn_segid==block1_cleaned.cellIDs(i) & df_cleaned.postsyn_segid == block1_cleaned.cellIDs(j));
            tempy = df_cleaned.postsyn_y(df_cleaned.presyn_segid==block1_cleaned.cellIDs(i) & df_cleaned.postsyn_segid == block1_cleaned.cellIDs(j));
            tempz = df_cleaned.postsyn_z(df_cleaned.presyn_segid==block1_cleaned.cellIDs(i) & df_cleaned.postsyn_segid == block1_cleaned.cellIDs(j));
            block1_cleaned.b1tob1postSynCoord  = [block1_cleaned.b1tob1postSynCoord;[tempx,tempy,tempz]];
            if ~isempty(tempx)
                if isExistReRoot(block1_cleaned.cellIDs(j))
                    temp2 = PathLengthToCoordinate([tempx,tempy,tempz],block1_cleaned.cellIDs(j));
                    block1_cleaned.b1tob1postSynPathLength = [block1_cleaned.b1tob1postSynPathLength;temp2];
                    normPathLength =  temp2 ./max(Pvec_tree(block1_cleaned.tree{j}));
                    block1_cleaned.b1tob1postSynPathLength_norm = [block1_cleaned.b1tob1postSynPathLength_norm;normPathLength];
                    clear normPathLength
                    clear temp2
                end
            end
            clear temp;
            clear tempx;
            clear tempy;
            clear tempz;
        end
    end
end


for i= 1:length(block2_cleaned.cellIDs)
    if isExistReRoot(block2_cleaned.cellIDs(i))
        block2_cleaned.tree(i) = SwctoZbrian(block2_cleaned.cellIDs(i));
    end
end

block1_cleaned.b1tob2size = [];
block1_cleaned.b1tob2postSynCoord = [];
block1_cleaned.b1tob2postSynPathLength = [];
block1_cleaned.b1tob2postSynPathLength_norm =[];

for i = 1:length(block1_cleaned.cellIDs) % axon
    for j = 1:length(block2_cleaned.cellIDs) % dendrite
        if ~isequal(block1_cleaned.cellIDs(i),block2_cleaned.cellIDs(j))
            temp = df_cleaned.size(df_cleaned.presyn_segid==block1_cleaned.cellIDs(i) & df_cleaned.postsyn_segid == block2_cleaned.cellIDs(j));
            block1_cleaned.b1tob2size  = [block1_cleaned.b1tob2size ;temp];
            clear temp;
            tempx = df_cleaned.postsyn_x(df_cleaned.presyn_segid==block1_cleaned.cellIDs(i) & df_cleaned.postsyn_segid == block2_cleaned.cellIDs(j));
            tempy = df_cleaned.postsyn_y(df_cleaned.presyn_segid==block1_cleaned.cellIDs(i) & df_cleaned.postsyn_segid == block2_cleaned.cellIDs(j));
            tempz = df_cleaned.postsyn_z(df_cleaned.presyn_segid==block1_cleaned.cellIDs(i) & df_cleaned.postsyn_segid == block2_cleaned.cellIDs(j));
            block1_cleaned.b1tob2postSynCoord  = [block1_cleaned.b1tob2postSynCoord;[tempx,tempy,tempz]];
            if ~isempty(tempx)
                if isExistReRoot(block2_cleaned.cellIDs(j))
                    temp2 = PathLengthToCoordinate([tempx,tempy,tempz],block2_cleaned.cellIDs(j));
                    block1_cleaned.b1tob2postSynPathLength = [block1_cleaned.b1tob2postSynPathLength; temp2 ];
                    normPathLength =  temp2./max(Pvec_tree(block2_cleaned.tree{j}));
                    block1_cleaned.b1tob2postSynPathLength_norm = [block1_cleaned.b1tob2postSynPathLength_norm;normPathLength];
                    clear normPathLength
                    clear temp2
                end
            end
            clear temp;
            clear tempx;
            clear tempy;
        end
    end
end

% mean synapse size in block2_cleaned_gamma038 (i-->j)
block2_cleaned.b2tob2size = [];
block2_cleaned.b2tob2postSynCoord = [];
block2_cleaned.b2tob2postSynPathLength = [];
block2_cleaned.b2tob2postSynPathLength_norm=[];

for i = 1:length(block2_cleaned.cellIDs) % axon
    for j = 1:length(block2_cleaned.cellIDs) % dendrite
        if ~isequal(block2_cleaned.cellIDs(i),block2_cleaned.cellIDs(j))
            temp = df_cleaned.size(df_cleaned.presyn_segid==block2_cleaned.cellIDs(i) & df_cleaned.postsyn_segid == block2_cleaned.cellIDs(j));
            block2_cleaned.b2tob2size = [block2_cleaned.b2tob2size;temp];
            clear temp;
            tempx = df_cleaned.postsyn_x(df_cleaned.presyn_segid==block2_cleaned.cellIDs(i) & df_cleaned.postsyn_segid == block2_cleaned.cellIDs(j));
            tempy = df_cleaned.postsyn_y(df_cleaned.presyn_segid==block2_cleaned.cellIDs(i) & df_cleaned.postsyn_segid == block2_cleaned.cellIDs(j));
            tempz = df_cleaned.postsyn_z(df_cleaned.presyn_segid==block2_cleaned.cellIDs(i) & df_cleaned.postsyn_segid == block2_cleaned.cellIDs(j));
            block2_cleaned.b2tob2postSynCoord  = [block2_cleaned.b2tob2postSynCoord;[tempx,tempy,tempz]];
            if ~isempty(tempx)
                if isExistReRoot(block2_cleaned.cellIDs(j))
                    temp2 = PathLengthToCoordinate([tempx,tempy,tempz],block2_cleaned.cellIDs(j));
                    block2_cleaned.b2tob2postSynPathLength = [block2_cleaned.b2tob2postSynPathLength;temp2];
                    normPathLength =   temp2./max(Pvec_tree(block2_cleaned.tree{j}));
                    block2_cleaned.b2tob2postSynPathLength_norm = [block2_cleaned.b2tob2postSynPathLength_norm;normPathLength];
                    clear normPathLength
                    clear temp2
                end
            end
            clear temp;
            clear tempx;
            clear tempy;
        end
    end
end


block2_cleaned.b2tob1size = [];
block2_cleaned.b2tob1postSynCoord = [];
block2_cleaned.b2tob1postSynPathLength =[];
block2_cleaned.b2tob1postSynPathLength_norm = [];


for i = 1:length(block2_cleaned.cellIDs) % axon
    for j = 1:length(block1_cleaned.cellIDs) % dendrite
        if ~isequal(block2_cleaned.cellIDs(i),block1_cleaned.cellIDs(j))
            temp = df_cleaned.size(df_cleaned.presyn_segid==block2_cleaned.cellIDs(i) & df_cleaned.postsyn_segid == block1_cleaned.cellIDs(j));
            block2_cleaned.b2tob1size = [block2_cleaned.b2tob1size;temp];
            clear temp;
            tempx = df_cleaned.postsyn_x(df_cleaned.presyn_segid==block2_cleaned.cellIDs(i) & df_cleaned.postsyn_segid == block1_cleaned.cellIDs(j));
            tempy = df_cleaned.postsyn_y(df_cleaned.presyn_segid==block2_cleaned.cellIDs(i) & df_cleaned.postsyn_segid == block1_cleaned.cellIDs(j));
            tempz = df_cleaned.postsyn_z(df_cleaned.presyn_segid==block2_cleaned.cellIDs(i) & df_cleaned.postsyn_segid == block1_cleaned.cellIDs(j));
            block2_cleaned.b2tob1postSynCoord  = [block2_cleaned.b2tob1postSynCoord;[tempx,tempy,tempz]];
            if ~isempty(tempx)
                if isExistReRoot(block1_cleaned.cellIDs(j))
                    temp2 = PathLengthToCoordinate([tempx,tempy,tempz],block1_cleaned.cellIDs(j));
                    block2_cleaned.b2tob1postSynPathLength = [block2_cleaned.b2tob1postSynPathLength;temp2];
                    normPathLength =   temp2./max(Pvec_tree(block1_cleaned.tree{j}));
                    block2_cleaned.b2tob1postSynPathLength_norm = [block2_cleaned.b2tob1postSynPathLength_norm;normPathLength];
                    clear normPathLength
                    clear temp2;
                end
            end
            clear tempz;
            clear tempx;
            clear tempy;
        end
    end
end

%%
figure;
subplot(4,4,1)
histogram(log10(block1_cleaned.b1tob1size),'BinWidth',0.05,'EdgeColor',[1,0.5,0],'DisplayStyle','stairs','LineWidth',2);
hold on;
histogram(log10(block2_cleaned.b2tob2size),'BinWidth',0.05,'EdgeColor',[0,0.5,1],'DisplayStyle','stairs','LineWidth',2);
%axis square;
box off;
set(gca,'XTick',[2,3,4],'XTickLabel',[10^2,10^3,10^4]);
xlabel('log size (vx)');
ylabel('count');
legend({'modA->modA','modO-modO'},'Location','bestoutside');
[h,p] = kstest2(log10(block1_cleaned.b1tob1size),log10(block2_cleaned.b2tob2size));
title(p);
% line([mean(log10(block1_cleaned.b1tob1size)),mean(log10(block1_cleaned.b1tob1size))],...
%     [0,500],'color',[1,0.5,0],'lineWidth',2);
% line([mean(log10(block2_cleaned_gamma038.b2tob2size)),mean(log10(block2_cleaned_gamma038.b2tob2size))],...
%     [0,500],'color',[0,0.5,1],'lineWidth',2);


subplot(4,4,2)
histogram(log10(block1_cleaned.b1tob2size),'BinWidth',0.05,'EdgeColor',[1,0.5,0],'DisplayStyle','stairs','LineWidth',2,'EdgeAlpha',0.4);
hold on;
histogram(log10(block2_cleaned.b2tob1size),'BinWidth',0.05,'EdgeColor',[0,0.5,1],'DisplayStyle','stairs','LineWidth',2,'EdgeAlpha',0.4);
%axis square;
box off;
set(gca,'XTick',[2,3,4],'XTickLabel',[10^2,10^3,10^4]);
legend({'modA->modO','modO->modA'},'Location','bestoutside');
[h,p] = kstest2(log10(block1_cleaned.b1tob2size),log10(block2_cleaned.b2tob1size));
title(p);

subplot(4,4,5)
histogram(block1_cleaned.b1tob1postSynPathLength_norm,'BinWidth',0.05,'EdgeColor','none','FaceColor',block1_cleaned.color,'Normalization','probability');
hold on
histogram(block2_cleaned.b2tob1postSynPathLength_norm,'BinWidth',0.05,'EdgeColor',block2_cleaned.color,'FaceColor','none','LineWidth',2,'DisplayStyle','stairs','Normalization','probability');
%axis square;
box off;
xlabel('Pathlength (\mum)');
ylabel('Relative probability');
legend({'modA->modA','modO->modA'},'Location','bestoutside');


subplot(4,4,6)
histogram(block2_cleaned.b2tob2postSynPathLength_norm,'BinWidth',0.05,'EdgeColor','none','FaceColor',block2_cleaned.color,'Normalization','probability');
hold on
histogram(block1_cleaned.b1tob2postSynPathLength_norm,'BinWidth',0.05,'EdgeColor',block1_cleaned.color,'FaceColor','none','LineWidth',2,'DisplayStyle','stairs','Normalization','probability');
%axis square;
box off;
xlabel('Pathlength (\mum)');
ylabel('Relative probability');
legend({'modO->modO','modA->modO'},'Location','bestoutside');



stats_mods_size = {'modA-->modA';'modA-->modO';'modO-->modO'; 'modO-->modA'};
means_mods_size = [mean(block1_cleaned.b1tob1size);mean(block1_cleaned.b1tob2size);...
         mean(block2_cleaned.b2tob2size);mean(block2_cleaned.b2tob1size)];
sdevs_mods_size =  [std(block1_cleaned.b1tob1size);std(block1_cleaned.b1tob2size);...
         std(block2_cleaned.b2tob2size);std(block2_cleaned.b2tob1size)];
nums_mods_size =  [numel(block1_cleaned.b1tob1size);numel(block1_cleaned.b1tob2size);...
          numel(block2_cleaned.b2tob2size);numel(block2_cleaned.b2tob1size)];


block2Table_size = table(means_mods_size,sdevs_mods_size,'RowNames',stats_mods_size);


stats_mods_pathlength = {'modA-->modA';'modA-->modO';'modO-->modO'; 'modO-->modA'};
means_mods_pathlength = [mean(block1_cleaned.b1tob1postSynPathLength);mean(block1_cleaned.b1tob2postSynPathLength);...
         mean(block2_cleaned.b2tob2postSynPathLength);mean(block2_cleaned.b2tob1postSynPathLength)];
sdevs_mods_pathlength =  [std(block1_cleaned.b1tob1postSynPathLength);std(block1_cleaned.b1tob2postSynPathLength);...
         std(block2_cleaned.b2tob2postSynPathLength);std(block2_cleaned.b2tob1postSynPathLength)];
nums_mods_pathlength =  [numel(block1_cleaned.b1tob1postSynPathLength);numel(block1_cleaned.b1tob2postSynPathLength);...
         numel(block2_cleaned.b2tob2postSynPathLength);numel(block2_cleaned.b2tob1postSynPathLength)];

block2Table_pathlengths = table(means_mods_pathlength,sdevs_mods_pathlength,'RowNames',stats_mods_pathlength);


subplot(4,4,3)
heatmap(block2Table_size{:,:},'YDisplayLabels',block2Table_size.Properties.RowNames,...
    'XDisplayLabels',block2Table_size.Properties.VariableNames,'FontSize',15)

% errorbar(means_mods_size,sdevs_mods_size./sqrt(nums_mods_size),'ko','lineWidth',2);
% box off;
% set(gca,'XTick',[1,2,3,4],'XTickLabels',stats_mods_size,'XTickLabelRotation',45);
% daspect([1,40,1]);


subplot(4,4,8)
heatmap(block2Table_pathlengths{:,:},'YDisplayLabels',block2Table_pathlengths.Properties.RowNames,...
    'XDisplayLabels',block2Table_pathlengths.Properties.VariableNames,'FontSize',15)
% 
% errorbar(means_mods_pathlength,sdevs_mods_pathlength./sqrt(nums_mods_pathlength),'ko','lineWidth',2);
% box off;
% set(gca,'XTick',[1,2,3,4],'XTickLabels',stats_mods_size,'XTickLabelRotation',45);
% daspect([1,4,1]);

%%

% mean synapse size in block2_cleaned_gamma038 to motor (i-->j)
ABDr_CellIDs = [82140,82145,82143,77648,82146,77302,77301,82194,77705,77710,77305,77672,77300,77709,77661];
ABDc_CellIDs = [82213,77654,77646,82212,77295,81172,77657,82197,77682,77154, 77652,77658,77628,77292,77688,82195,77296];
ABDIr_CellIDs = [78553,77150,77886,78547,77631,77158,77665,78552,77668,77618,77634,78556];
ABDIc_CellIDs = [78574,79051,77148,77625,77144,77692,77641,79066,77640,77643];

ABDmotor = [ABDr_CellIDs,ABDc_CellIDs];
ABDinter  = [ABDIr_CellIDs,ABDIc_CellIDs];


for i= 1:length(block1_cleaned_submod.cellIDs)
    if isExistReRoot(block1_cleaned_submod.cellIDs(i))
        block1_cleaned_submod.tree(i) = SwctoZbrian(block1_cleaned_submod.cellIDs(i));
    end
end

for i= 1:length(block2_cleaned_submod.cellIDs)
    if isExistReRoot(block2_cleaned_submod.cellIDs(i))
        block2_cleaned_submod.tree(i) = SwctoZbrian(block2_cleaned_submod.cellIDs(i));
    end
end

block1_cleaned_submod.b1tob1 = [];
block1_cleaned_submod.b1tob1PathLength = [];
block1_cleaned_submod.b1tob1PostSynCoord = [];
block1_cleaned_submod.b1tob1PathLength_norm =[];

for i = 1:length(block1_cleaned_submod.cellIDs) % axon
    for j = 1:length(block1_cleaned_submod.cellIDs) % dendrite
        if ~isequal(block1_cleaned_submod.cellIDs(i),block1_cleaned_submod.cellIDs(j))
            temp = df_cleaned.size(df_cleaned.presyn_segid==block1_cleaned_submod.cellIDs(i) & df_cleaned.postsyn_segid == block1_cleaned_submod.cellIDs(j));
            block1_cleaned_submod.b1tob1 = [block1_cleaned_submod.b1tob1;temp];
            clear temp;
            tempx = df_cleaned.postsyn_x(df_cleaned.presyn_segid==block1_cleaned_submod.cellIDs(i) & df_cleaned.postsyn_segid == block1_cleaned_submod.cellIDs(j));
            tempy = df_cleaned.postsyn_y(df_cleaned.presyn_segid==block1_cleaned_submod.cellIDs(i) & df_cleaned.postsyn_segid == block1_cleaned_submod.cellIDs(j));
            tempz = df_cleaned.postsyn_z(df_cleaned.presyn_segid==block1_cleaned_submod.cellIDs(i) & df_cleaned.postsyn_segid == block1_cleaned_submod.cellIDs(j));
            block1_cleaned_submod.b1tob1PostSynCoord  = [block1_cleaned_submod.b1tob1PostSynCoord;[tempx,tempy,tempz]];
            if ~isempty(tempx)
                if isExistReRoot(block1_cleaned_submod.cellIDs(j))
                    block1_cleaned_submod.b1tob1PathLength = [block1_cleaned_submod.b1tob1PathLength;PathLengthToCoordinate([tempx,tempy,tempz],block1_cleaned_submod.cellIDs(j))];
                     normPathLength =   PathLengthToCoordinate([tempx,tempy,tempz],block1_cleaned_submod.cellIDs(j))./max(Pvec_tree(block1_cleaned_submod.tree{j}));
                     block1_cleaned_submod.b1tob1PathLength_norm = [block1_cleaned_submod.b1tob1PathLength_norm; normPathLength];
                end
            end
            clear tempz;
            clear tempx;
            clear tempy;
            
        end
    end
end

block1_cleaned_submod.b1tob2 = [];
block1_cleaned_submod.b1tob2PathLength = [];
block1_cleaned_submod.b1tob2PostSynCoord = [];
 block1_cleaned_submod.b1tob2PathLength_norm = [];

for i = 1:length(block1_cleaned_submod.cellIDs) % axon
    for j = 1:length(block2_cleaned_submod.cellIDs) % dendrite
        if ~isequal(block1_cleaned_submod.cellIDs(i),block2_cleaned_submod.cellIDs(j))
            temp = df_cleaned.size(df_cleaned.presyn_segid==block1_cleaned_submod.cellIDs(i) & df_cleaned.postsyn_segid == block2_cleaned_submod.cellIDs(j));
            block1_cleaned_submod.b1tob2 = [block1_cleaned_submod.b1tob2;temp];
            clear temp;
            tempx = df_cleaned.postsyn_x(df_cleaned.presyn_segid==block1_cleaned_submod.cellIDs(i) & df_cleaned.postsyn_segid == block2_cleaned_submod.cellIDs(j));
            tempy = df_cleaned.postsyn_y(df_cleaned.presyn_segid==block1_cleaned_submod.cellIDs(i) & df_cleaned.postsyn_segid == block2_cleaned_submod.cellIDs(j));
            tempz = df_cleaned.postsyn_z(df_cleaned.presyn_segid==block1_cleaned_submod.cellIDs(i) & df_cleaned.postsyn_segid == block2_cleaned_submod.cellIDs(j));
            block1_cleaned_submod.b1tob2PostSynCoord  = [block1_cleaned_submod.b1tob2PostSynCoord;[tempx,tempy,tempz]];
            if ~isempty(tempx)
                if isExistReRoot(block2_cleaned_submod.cellIDs(j))
                     block1_cleaned_submod.b1tob2PathLength = [block1_cleaned_submod.b1tob2PathLength;PathLengthToCoordinate([tempx,tempy,tempz],block1_cleaned_submod.cellIDs(j))];
                     normPathLength =   PathLengthToCoordinate([tempx,tempy,tempz],block1_cleaned_submod.cellIDs(j))./max(Pvec_tree(block2_cleaned_submod.tree{j}));
                     block1_cleaned_submod.b1tob2PathLength_norm = [block1_cleaned_submod.b1tob2PathLength_norm ; normPathLength];
                end
            end
            clear tempz;
            clear tempx;
            clear tempy;
            
        end
    end
end



block2_cleaned_submod.b2tob2 = [];
block2_cleaned_submod.b2tob2PathLength = [];
block2_cleaned_submod.b2tob2PostSynCoord = [];
block2_cleaned_submod.b2tob2PathLength_norm = [];



for i = 1:length(block2_cleaned_submod.cellIDs) % axon
    for j = 1:length(block2_cleaned_submod.cellIDs) % dendrite
        if ~isequal(block2_cleaned_submod.cellIDs(i),block2_cleaned_submod.cellIDs(j))
            temp = df_cleaned.size(df_cleaned.presyn_segid==block2_cleaned_submod.cellIDs(i) & df_cleaned.postsyn_segid == block2_cleaned_submod.cellIDs(j));
            block2_cleaned_submod.b2tob2 = [block2_cleaned_submod.b2tob2;temp];
            clear temp;
            tempx = df_cleaned.postsyn_x(df_cleaned.presyn_segid==block2_cleaned_submod.cellIDs(i) & df_cleaned.postsyn_segid == block2_cleaned_submod.cellIDs(j));
            tempy = df_cleaned.postsyn_y(df_cleaned.presyn_segid==block2_cleaned_submod.cellIDs(i) & df_cleaned.postsyn_segid == block2_cleaned_submod.cellIDs(j));
            tempz = df_cleaned.postsyn_z(df_cleaned.presyn_segid==block2_cleaned_submod.cellIDs(i) & df_cleaned.postsyn_segid == block2_cleaned_submod.cellIDs(j));
            block2_cleaned_submod.b2tob2PostSynCoord  = [block2_cleaned_submod.b2tob2PostSynCoord;[tempx,tempy,tempz]];
            if ~isempty(tempx)
                if isExistReRoot(block2_cleaned_submod.cellIDs(j))
                     block2_cleaned_submod.b2tob2PathLength = [block2_cleaned_submod.b2tob2PathLength;PathLengthToCoordinate([tempx,tempy,tempz],block2_cleaned_submod.cellIDs(j))];
                     normPathLength =   PathLengthToCoordinate([tempx,tempy,tempz],block2_cleaned_submod.cellIDs(j))./max(Pvec_tree(block2_cleaned_submod.tree{j}));
                     block2_cleaned_submod.b2tob2PathLength_norm = [ block2_cleaned_submod.b2tob2PathLength_norm ; normPathLength];
                end
            end
            clear tempz;
            clear tempx;
            clear tempy;
            
        end
    end
end


block2_cleaned_submod.b2tob1 = [];
block2_cleaned_submod.b2tob1PathLength = [];
block2_cleaned_submod.b2tob1PostSynCoord = [];
block2_cleaned_submod.b2tob1PathLength_norm = [];

for i = 1:length(block2_cleaned_submod.cellIDs) % axon
    for j = 1:length(block1_cleaned_submod.cellIDs) % dendrite
        if ~isequal(block2_cleaned_submod.cellIDs(i),block1_cleaned_submod.cellIDs(j))
            temp = df_cleaned.size(df_cleaned.presyn_segid==block2_cleaned_submod.cellIDs(i) & df_cleaned.postsyn_segid == block1_cleaned_submod.cellIDs(j));
            block2_cleaned_submod.b2tob1 = [block2_cleaned_submod.b2tob1;temp];
            clear temp;
            tempx = df_cleaned.postsyn_x(df_cleaned.presyn_segid==block2_cleaned_submod.cellIDs(i) & df_cleaned.postsyn_segid == block1_cleaned_submod.cellIDs(j));
            tempy = df_cleaned.postsyn_y(df_cleaned.presyn_segid==block2_cleaned_submod.cellIDs(i) & df_cleaned.postsyn_segid == block1_cleaned_submod.cellIDs(j));
            tempz = df_cleaned.postsyn_z(df_cleaned.presyn_segid==block2_cleaned_submod.cellIDs(i) & df_cleaned.postsyn_segid == block1_cleaned_submod.cellIDs(j));
            block2_cleaned_submod.b2tob1PostSynCoord  = [block2_cleaned_submod.b2tob1PostSynCoord;[tempx,tempy,tempz]];
            if ~isempty(tempx)
                if isExistReRoot(block1_cleaned_submod.cellIDs(j))
                     block2_cleaned_submod.b2tob1PathLength = [block2_cleaned_submod.b2tob1PathLength;PathLengthToCoordinate([tempx,tempy,tempz],block1_cleaned_submod.cellIDs(j))];
                     normPathLength =    PathLengthToCoordinate([tempx,tempy,tempz],block1_cleaned_submod.cellIDs(j))./max(Pvec_tree(block1_cleaned_submod.tree{j}));
                     block2_cleaned_submod.b2tob1PathLength_norm = [block2_cleaned_submod.b2tob1PathLength_norm ; normPathLength];
                end
            end
            clear tempz;
            clear tempx;
            clear tempy;
            
        end
    end
end

%%

figure;
subplot(4,4,1)
histogram(log10(block1_cleaned_submod.b1tob1),'BinWidth',0.05,'EdgeColor',block1_cleaned_submod.color,'DisplayStyle','stairs','LineWidth',2);
hold on;
histogram(log10(block2_cleaned_submod.b2tob2),'BinWidth',0.05,'EdgeColor',block2_cleaned_submod.color,'DisplayStyle','stairs','LineWidth',2);
%axis square;
box off;
set(gca,'XTick',[2,3,4],'XTickLabel',[10^2,10^3,10^4]);
xlabel('log size (vx)');
ylabel('count');
legend({'modO_I->modO_I','modO_M-modO_M'},'Location','bestoutside');
[h,p] = kstest2(log10(block1_cleaned_submod.b1tob1),log10(block2_cleaned_submod.b2tob2));
title(p);


subplot(4,4,2)
histogram(log10(block1_cleaned_submod.b1tob2),'BinWidth',0.05,'EdgeColor',block1_cleaned_submod.color,'DisplayStyle','stairs','LineWidth',2,'EdgeAlpha',0.4);
hold on;
histogram(log10(block2_cleaned_submod.b2tob1),'BinWidth',0.05,'EdgeColor',block2_cleaned_submod.color,'DisplayStyle','stairs','LineWidth',2,'EdgeAlpha',0.4);
%axis square;
box off;
set(gca,'XTick',[2,3,4],'XTickLabel',[10^2,10^3,10^4]);
legend({'modO_I->modO_M','modO_M->modO_I'},'Location','bestoutside');
[h,p] = kstest2(log10(block1_cleaned_submod.b1tob2),log10(block2_cleaned_submod.b2tob1));
title(p);


subplot(4,4,5)
histogram(block1_cleaned_submod.b1tob1PathLength,'BinWidth',20,'EdgeColor','none','FaceColor',block1_cleaned_submod.color,'Normalization','probability');
hold on
histogram(block2_cleaned_submod.b2tob1PathLength,'BinWidth',20,'EdgeColor',block2_cleaned_submod.color,'FaceColor','none','LineWidth',2,'DisplayStyle','stairs','Normalization','probability');
%axis square;
box off;
xlabel('Pathlength (\mum)');
ylabel('count');
legend({'modO_I->modO_I','modO_M->modO_I'},'Location','bestoutside');


subplot(4,4,6)
histogram(block2_cleaned_submod.b2tob2PathLength,'BinWidth',20,'EdgeColor','none','FaceColor',block2_cleaned_submod.color,'Normalization','probability');
hold on
histogram(block1_cleaned_submod.b1tob2PathLength,'BinWidth',20,'EdgeColor',block1_cleaned_submod.color,'FaceColor','none','LineWidth',2,'DisplayStyle','stairs','Normalization','probability');
%axis square;
box off;
xlabel('Pathlength (\mum)');
ylabel('count');
legend({'modO_M->modO_M','modO_I->modO_M'},'Location','bestoutside');


%subplot(4,4,3)

stats_submods_pathlengths = {'modOi-->modOi';'modOi-->modOm';'modOm-->modOm'; 'modOm-->modOi'};
means_submods_pathlengths = [mean(block1_cleaned_submod.b1tob1PathLength);mean(block1_cleaned_submod.b1tob2PathLength);...
         mean(block2_cleaned_submod.b2tob2PathLength);mean(block2_cleaned_submod.b2tob1PathLength)];
sdevs_submods_pathlengths =  [std(block1_cleaned_submod.b1tob1PathLength);std(block1_cleaned_submod.b1tob2PathLength);...
         std(block2_cleaned_submod.b2tob2PathLength);std(block2_cleaned_submod.b2tob1PathLength)];

block2Table_submods_pathlengths = table(means_submods_pathlengths,sdevs_submods_pathlengths,'RowNames',stats_submods_pathlengths);



stats_submods_size = {'modOi-->modOi';'modOi-->modOm';'modOm-->modOm'; 'modOm-->modOi'};
means_submods_size = [mean(block1_cleaned_submod.b1tob1);mean(block1_cleaned_submod.b1tob2);...
         mean(block2_cleaned_submod.b2tob2);mean(block2_cleaned_submod.b2tob1)];
sdevs_submods_size =  [std(block1_cleaned_submod.b1tob1);std(block1_cleaned_submod.b1tob2);...
         std(block2_cleaned_submod.b2tob2);std(block2_cleaned_submod.b2tob1)];

block2Table_submods_size = table(means_submods_size,sdevs_submods_size,'RowNames',stats_submods_size);


subplot(4,4,3)
heatmap(block2Table_submods_size{:,:},'YDisplayLabels',block2Table_submods_size.Properties.RowNames,...
    'XDisplayLabels',block2Table_submods_size.Properties.VariableNames)


subplot(4,4,7)
heatmap(block2Table_submods_pathlengths{:,:},'YDisplayLabels',block2Table_submods_pathlengths.Properties.RowNames,...
    'XDisplayLabels',block2Table_submods_pathlengths.Properties.VariableNames)
%%

for i = 1:length(block1_cleaned.cellIDs)  
    if ~isempty(block1_cleaned.tree{i})
   block1_cleaned.treeLength(i) =  sum(len_tree(block1_cleaned.tree{i}));
    end
end

block1_cleaned.treeLength(block1_cleaned.treeLength==0) = [];

for i = 1:length(block2_cleaned.cellIDs)  
    if ~isempty(block2_cleaned.tree{i})
   block2_cleaned.treeLength(i) =  sum(len_tree(block2_cleaned.tree{i}));
    end
end

 block2_cleaned.treeLength( block2_cleaned.treeLength==0) = [];

 figure;
 subplot(4,4,6)
 histogram(rmoutliers(block1_cleaned.treeLength),'BinWidth',80,'EdgeColor',block1_cleaned.color,'LineWidth',2,'DisplayStyle','stairs')
 hold on
 histogram(rmoutliers(block2_cleaned.treeLength),'BinWidth',80,'EdgeColor',block2_cleaned.color,'LineWidth',2,'DisplayStyle','stairs')
 xlabel('Pathlength (\mum)');
ylabel('frequency');
legend({'modA','modO'},'Location','bestoutside');
box off;

%%
figure;

subplot(4,4,1)
histogram(log10(block2_cleaned.b2toABDm),'BinWidth',0.1,'EdgeColor',block2_cleaned_submod.color,'DisplayStyle','stairs','LineWidth',2);
hold on;
histogram(log10(block2_cleaned.b2toABDi),'BinWidth',0.1,'EdgeColor',block1_cleaned_submod.color,'DisplayStyle','stairs','LineWidth',2);
box off;
axis square;
set(gca,'XTick',[2,3,4],'XTickLabel',[10^2,10^3,10^4]);
legend({'mod2b->ABD_m','mod2a-> ABD_i'},'Location','bestoutside');
xlabel('log size (vx)');
ylabel('count');
[h,p] = kstest2(log10(block2_cleaned.b2toABDm),log10(block2_cleaned.b2toABDi));
title(p);

subplot(4,4,2)
histogram(log10(block2_cleaned.b2tob2size),'BinWidth',0.1,'EdgeColor',[0,0.5,1],'DisplayStyle','stairs','LineWidth',2);
hold on;
histogram(log10(block2_cleaned.b2toABDm),'BinWidth',0.1,'EdgeColor','g','DisplayStyle','stairs','LineWidth',2);
histogram(log10(block2_cleaned.b2toABDi),'BinWidth',0.1,'EdgeColor','m','DisplayStyle','stairs','LineWidth',2);
box off;
axis square;
set(gca,'XTick',[2,3,4],'XTickLabel',[10^2,10^3,10^4]);
legend({'b2->b2','b2->ABD_m','b2-> ABD_i'},'Location','bestoutside');


%%
figure;
subplot(4,4,1)
histogram(block1_cleaned.b1tob1postSynPathLength,'BinWidth',20,'EdgeColor','none','FaceColor',[1,0.5,0],'Normalization','probability');
hold on
histogram(block2_cleaned.b2tob1postSynPathLength,'BinWidth',20,'EdgeColor',[0,0.5,1],'FaceColor','none','LineWidth',2,'DisplayStyle','stairs','Normalization','probability');

%axis square;
box off;
xlabel('Pathlength (\mum)');
ylabel('count');
legend({'mod1->mod1','mod2->mod1'},'Location','bestoutside');


subplot(4,4,3)
histogram(block2_cleaned.b2tob2postSynPathLength,'BinWidth',20,'EdgeColor','none','FaceColor',[0,0.5,1],'Normalization','probability');
hold on;
histogram(block1_cleaned.b1tob2postSynPathLength,'BinWidth',20,'EdgeColor',[1,0.5,0],'FaceColor','none','LineWidth',2,'DisplayStyle','stairs','Normalization','probability');

%axis square;
box off;
xlabel('Pathlength (\mum)');
ylabel('count');
legend({'mod2->mod2','mod1->mod2'},'Location','bestoutside');

subplot(4,4,5)
histogram(block2_cleaned.b2toABDmPathLength,'BinWidth',20,'EdgeColor','none','EdgeColor',block2_cleaned_submod.color,...
    'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
hold on
histogram(block2_cleaned.b2toABDiPathLength,'BinWidth',20,'EdgeColor','none','EdgeColor',block1_cleaned_submod.color,...
    'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
axis square;
box off;
xlabel('Pathlength (\mum)');
ylabel('count');
[h,p] = kstest2(block2_cleaned.b2toABDmPathLength,block2_cleaned.b2toABDiPathLength);
title(p);

%% plotSpreads

withinModule = [rmoutliers(block1_cleaned.b1tob1postSynPathLength);...
    rmoutliers(block2_cleaned.b2tob1postSynPathLength)];
ids = [ones(length(rmoutliers(block1_cleaned.b1tob1postSynPathLength)),1);...
    2*ones(length(rmoutliers(block2_cleaned.b2tob1postSynPathLength)),1)];

g = gramm('x',ids,'y',withinModule);
%g.geom_jitter();
g.stat_violin()
g.draw;

%% dual synapses coorelations
block1_cleaned.pairs_all = [];
block1_cleaned.pairs_single = [];
block1_cleaned.pairs_dual = [];
block1_cleaned.pairs_triple = [];
block1_cleaned.pairs_quart = [];
block1_cleaned.pairs_pent = [];

for i = 1:length(block1_cleaned.cellIDs)
    for j = 1:length(block1_cleaned.cellIDs)
        if ~isequal(block1_cleaned.cellIDs(i),block1_cleaned.cellIDs(j))
            sz = df_cleaned.size(df_cleaned.presyn_segid==block1_cleaned.cellIDs(i) & df_cleaned.postsyn_segid == block1_cleaned.cellIDs(j));
            block1_cleaned.pairs_all = [block1_cleaned.pairs_all;repmat(block1_cleaned.cellIDs(i),length(sz),1),repmat(block1_cleaned.cellIDs(j),length(sz),1),sz];
            if size(sz,1) == 1
                block1_cleaned.pairs_single =  [block1_cleaned.pairs_single;repmat(block1_cleaned.cellIDs(i),length(sz),1),repmat(block1_cleaned.cellIDs(j),length(sz),1),sz];
            elseif size(sz,1) == 2
                block1_cleaned.pairs_dual =  [block1_cleaned.pairs_dual;repmat(block1_cleaned.cellIDs(i),length(sz),1),repmat(block1_cleaned.cellIDs(j),length(sz),1),sz];
            elseif size(sz,1) == 3
                block1_cleaned.pairs_triple =  [block1_cleaned.pairs_triple;repmat(block1_cleaned.cellIDs(i),length(sz),1),repmat(block1_cleaned.cellIDs(j),length(sz),1),sz];
            elseif size(sz,1) == 4
                block1_cleaned.pairs_quart =  [block1_cleaned.pairs_quart;repmat(block1_cleaned.cellIDs(i),length(sz),1),repmat(block1_cleaned.cellIDs(j),length(sz),1),sz];
            else
                block1_cleaned.pairs_pent =  [block1_cleaned.pairs_pent;repmat(block1_cleaned.cellIDs(i),length(sz),1),repmat(block1_cleaned.cellIDs(j),length(sz),1),sz];
            end
            clear sz
        end
    end
end

%scatter(block1_cleaned.pairs_dual(1:2:end,3),block1_cleaned.pairs_dual(2:2:end,3))
block1_cleaned.dualpairs = [];
for i = 1:length(block1_cleaned.pairs_dual)
    block1_cleaned.dualpairs = [block1_cleaned.dualpairs;nchoosek(block1_cleaned.pairs_dual(i:i+1,3),2)];
end
block1_cleaned.dualtriple = [];
for i = 1:length(block1_cleaned.pairs_triple)
    block1_cleaned.dualtriple = [block1_cleaned.dualtriple;nchoosek(block1_cleaned.pairs_triple(i:i+2,3),2)];
end
block1_cleaned.dualquart = [];
for i = 1:length(block1_cleaned.pairs_quart)
    block1_cleaned.dualquart = [block1_cleaned.dualquart;nchoosek(block1_cleaned.pairs_quart(i:i+3,3),2)];
end
block1_cleaned.dualpent = [];
for i = 1:length(block1_cleaned.pairs_pent)
    block1_cleaned.dualpent = [block1_cleaned.dualpent;nchoosek(block1_cleaned.pairs_pent(i:i+4,3),2)];
end
% 
% block1_cleaned.diffPairs = abs(block1_cleaned.pairs_all(:,1)-block1_cleaned.pairs_all(:,2))./block1_cleaned.pairs_all(:,1); % all pairs
% block1_cleaned.DualdiffPairs = abs(block1_cleaned.pairs_dual(:,1)-block1_cleaned.pairs_dual(:,2)); % only dual pairs
% block1_cleaned.locs = diff(block1_cleaned.DualdiffPairs);
% block1_cleaned.index = find(block1_cleaned.locs);

% ind = 1;
% block1_cleaned.dualPairs = [];
% for i = 1:length(block1_cleaned.index)
%     nums = block1_cleaned.index(i)-ind;
%     block1_cleaned.dualPairs = [block1_cleaned.DualdiffPairs;nchoosek(block1_cleaned.pairs_dual(ind:block1_cleaned.index(i),3),2)];
%     ind = block1_cleaned.index(i)+1;
% end


% block1_cleaned-->b2
block1_cleaned.b1tob2pairs_all = [];
block1_cleaned.b1tob2pairs_dual = [];

for i = 1:length(block1_cleaned.cellIDs)
    for j = 1:length(block2_cleaned.cellIDs)
        if ~isequal(block1_cleaned.cellIDs(i),block2_cleaned.cellIDs(j))
            sz = df_cleaned.size(df_cleaned.presyn_segid==block1_cleaned.cellIDs(i) & df_cleaned.postsyn_segid == block2_cleaned.cellIDs(j));
            block1_cleaned.b1tob2pairs_all = [block1_cleaned.b1tob2pairs_all;repmat(block1_cleaned.cellIDs(i),length(sz),1),repmat(block2_cleaned.cellIDs(j),length(sz),1),sz];
            if size(sz,1)>1
                block1_cleaned.b1tob2pairs_dual =  [block1_cleaned.b1tob2pairs_dual;repmat(block1_cleaned.cellIDs(i),length(sz),1),repmat(block2_cleaned.cellIDs(j),length(sz),1),sz];
            end
            clear sz
        end
    end
end

block1_cleaned.b1tob2diffPairs = abs(block1_cleaned.b1tob2pairs_all(:,1)-block1_cleaned.b1tob2pairs_all(:,2))./block1_cleaned.b1tob2pairs_all(:,1);



% % gmm for dual pairs
% 
% % monovariate
% block1_cleaned.gmmDualpairs = fitgmdist(log10(block1_cleaned.pairs_dual(:,3)),2,'Options',statset('MaxIter',1000));
% samples = 2:0.01:5;
% evalat = samples(randi(length(samples),1000,1))';
% block1_cleaned.gmmDualpairsPdf_cleaned = pdf_cleaned(block1_cleaned.gmmDualpairs,[2:0.01:5]');
% 
% %bivariate
% 
% block1_cleaned.gmmDualpairsBI = fitgmdist([log10(block1_cleaned.dualPairs(:,1)),log10(block1_cleaned.dualPairs(:,2))],2,'Options',statset('MaxIter',1000));
% block1_cleaned.gmmDualpairsBIPdf_cleaned = @(x,y)reshape(pdf_cleaned(block1_cleaned.gmmDualpairsBI,[x(:),y(:)]),size(x));
% 
% 
% subplot(2,2,1);
% histogram(log10(block1_cleaned.pairs_dual(:,3)),'BinWidth',0.1,'EdgeColor',[1,0.5,0],'DisplayStyle','stairs','LineWidth',2);
% xlabel('size(vx)');
% ylabel('count');
% set(gca,'XTick',[2,3,4,5],'XTickLabels',[10^2,10^3,10^4,10^5]);
% hold on
% yyaxis right;
% plot(2:0.01:5,block1_cleaned.gmmDualpairsPdf_cleaned,'color',[1,0.5,0],'LineWidth',2);
% axis square;
% 
% subplot(2,2,2)
% scatter(log10(block1_cleaned.dualPairs(:,1)),log10(block1_cleaned.dualPairs(:,2)),'.','markerEdgecolor',[1,0.5,0],'MarkerEdgeAlpha',0.2);
% hold on;
% fcontour(block1_cleaned.gmmDualpairsBIPdf_cleaned,[2 5],'LineWidth',2);
% set(gca,'XTick',[2,3,4,5],'XTickLabels',[10^2,10^3,10^4,10^5],'YTick',[2,3,4,5],'YTickLabels',[10^2,10^3,10^4,10^5]);
% 
% axis square;

block2_cleaned.pairs_all = [];
block2_cleaned.pairs_single = [];
block2_cleaned.pairs_dual = [];
block2_cleaned.pairs_triple = [];
block2_cleaned.pairs_quart = [];
block2_cleaned.pairs_pent = [];

for i = 1:length(block2_cleaned.cellIDs)
    for j = 1:length(block2_cleaned.cellIDs)
        if ~isequal(block2_cleaned.cellIDs(i),block2_cleaned.cellIDs(j))
            sz = df_cleaned.size(df_cleaned.presyn_segid==block2_cleaned.cellIDs(i) & df_cleaned.postsyn_segid == block2_cleaned.cellIDs(j));
            block2_cleaned.pairs_all = [block2_cleaned.pairs_all;repmat(block2_cleaned.cellIDs(i),length(sz),1),repmat(block2_cleaned.cellIDs(j),length(sz),1),sz];
            if size(sz,1) == 1
                block2_cleaned.pairs_single =  [block2_cleaned.pairs_single;repmat(block2_cleaned.cellIDs(i),length(sz),1),repmat(block2_cleaned.cellIDs(j),length(sz),1),sz];
            elseif size(sz,1) == 2
                block2_cleaned.pairs_dual =  [block2_cleaned.pairs_dual;repmat(block2_cleaned.cellIDs(i),length(sz),1),repmat(block2_cleaned.cellIDs(j),length(sz),1),sz];
            elseif size(sz,1) == 3
                block2_cleaned.pairs_triple =  [block2_cleaned.pairs_triple;repmat(block2_cleaned.cellIDs(i),length(sz),1),repmat(block2_cleaned.cellIDs(j),length(sz),1),sz];
            elseif size(sz,1) == 4
                block2_cleaned.pairs_quart =  [block2_cleaned.pairs_quart;repmat(block2_cleaned.cellIDs(i),length(sz),1),repmat(block2_cleaned.cellIDs(j),length(sz),1),sz];
            else
                block2_cleaned.pairs_pent =  [block2_cleaned.pairs_pent;repmat(block2_cleaned.cellIDs(i),length(sz),1),repmat(block2_cleaned.cellIDs(j),length(sz),1),sz];
            end
            clear sz
        end
    end
end

block2_cleaned.dualpairs = [];
for i = 1:length(block1_cleaned.pairs_dual)
    block2_cleaned.dualpairs = [block2_cleaned.dualpairs;nchoosek(block2_cleaned.pairs_dual(i:i+1,3),2)];
end
block2_cleaned.dualtriple = [];
for i = 1:length(block2_cleaned.pairs_triple)
    block2_cleaned.dualtriple = [block2_cleaned.dualtriple;nchoosek(block2_cleaned.pairs_triple(i:i+2,3),2)];
end
block2_cleaned.dualquart = [];
for i = 1:length(block2_cleaned.pairs_quart)
    block2_cleaned.dualquart = [block2_cleaned.dualquart;nchoosek(block2_cleaned.pairs_quart(i:i+3,3),2)];
end
block2_cleaned.dualpent = [];
for i = 1:length(block2_cleaned.pairs_pent)
    block2_cleaned.dualpent = [block2_cleaned.dualpent;nchoosek(block2_cleaned.pairs_pent(i:i+4,3),2)];
end


% block2_cleaned_gamma038.diffPairs = abs(block2_cleaned_gamma038.pairs_all(:,1)-block2_cleaned_gamma038.pairs_all(:,2));
% block2_cleaned_gamma038.DualdiffPairs = abs(block2_cleaned_gamma038.pairs_dual(:,1)-block2_cleaned_gamma038.pairs_dual(:,2));
% 
% block2_cleaned_gamma038.locs = diff(block2_cleaned_gamma038.diffPairs);
% block2_cleaned_gamma038.index = find(block2_cleaned_gamma038.locs);
% 
% ind = 1;
% block2_cleaned_gamma038.dualPairs = [];
% for i = 1:length(block2_cleaned_gamma038.index)
%     nums = block2_cleaned_gamma038.index(i)-ind;
%     block2_cleaned_gamma038.dualPairs = [block2_cleaned_gamma038.dualPairs;nchoosek(block2_cleaned_gamma038.pairs_dual(ind:block2_cleaned_gamma038.index(i),3),2)];
%     ind = block2_cleaned_gamma038.index(i)+1;
% end
% 
%b2-->b1
block2_cleaned.b2tob1pairs_all = [];
block2_cleaned.b2tob1pairs_dual = [];

for i = 1:length(block2_cleaned.cellIDs)
    for j = 1:length(block1_cleaned.cellIDs)
        if ~isequal(block2_cleaned.cellIDs(i),block1_cleaned.cellIDs(j))
            sz = df_cleaned.size(df_cleaned.presyn_segid==block2_cleaned.cellIDs(i) & df_cleaned.postsyn_segid == block1_cleaned.cellIDs(j));
            block2_cleaned.b2tob1pairs_all = [block2_cleaned.b2tob1pairs_all;repmat(block2_cleaned.cellIDs(i),length(sz),1),repmat(block1_cleaned.cellIDs(j),length(sz),1),sz];
            if size(sz,1)>1
                block2_cleaned.b2tob1pairs_dual =  [block2_cleaned.b2tob1pairs_dual;repmat(block2_cleaned.cellIDs(i),length(sz),1),repmat(block1_cleaned.cellIDs(j),length(sz),1),sz];
            end
            clear sz
        end
    end
end

block2_cleaned.b2tob1diffPairs = abs(block2_cleaned.b2tob1pairs_all(:,1)-block2_cleaned.b2tob1pairs_all(:,2));

%block2_cleaned_gamma038.DualdiffPairs = abs(block2_cleaned_gamma038.pairs_dual(:,1)-block2_cleaned_gamma038.pairs_dual(:,2));



figure;
subplot(4,4,1)
scatter(log10(block1_cleaned.dualPairs(:,1)),log10(block1_cleaned.dualPairs(:,2)),'.','MarkerEdgeColor',[1,0.5,0],'MarkerEdgeAlpha',0.2);
hold on
f = showfit(ezfit(log10(block1_cleaned.dualPairs(:,1)),log10(block1_cleaned.dualPairs(:,2)),'affine'),'fitcolor',[1,0.5,0],'dispeqboxmode','off');
text(2.5,4,sprintf('r = %0.4f',f.r));
axis square;
set(gca,'XTick',[2,3,4,5],'XTickLabels',[10^2,10^3,10^4,10^5],'YTick',[2,3,4,5],'YTickLabels',[10^2,10^3,10^4,10^5]);

subplot(4,4,2)
scatter(log10(block2_cleaned.dualPairs(:,1)),log10(block2_cleaned.dualPairs(:,2)),'.','MarkerEdgeColor',[0,0.5,1],'MarkerEdgeAlpha',0.2);
hold on
f = showfit(ezfit(log10(block2_cleaned.dualPairs(:,1)),log10(block2_cleaned.dualPairs(:,2)),'affine'),'fitcolor',[0,0.5,1],'dispeqboxmode','off');
text(2.5,4,sprintf('r = %0.4f',f.r));
axis square;
set(gca,'XTick',[2,3,4,5],'XTickLabels',[10^2,10^3,10^4,10^5],'YTick',[2,3,4,5],'YTickLabels',[10^2,10^3,10^4,10^5]);

% h = scatterhist([log10(block1_cleaned.dualPairs(:,1));log10(block2_cleaned_gamma038.dualPairs(:,1))],...
%     [log10(block1_cleaned.dualPairs(:,2));log10(block2_cleaned_gamma038.dualPairs(:,2))],...
%     'Group',[ones(size(block1_cleaned.dualPairs,1),1);2*ones(size(block2_cleaned_gamma038.dualPairs,1),1)],'Kernel','on','color',[[1,0.5,0];[0,0.5,1]],...
%     'marker','.');
% set(h(1), 'XTick',[2,3,4,5],'XTickLabel',[10^2,10^3,10^4],'YTick',[2,3,4,5],'YTickLabel',[10^2,10^3,10^4]);
figure;
subplot(4,4,1)
histogram(histcounts(block1_cleaned.b1tob2diffPairs,unique( block1_cleaned.b1tob2diffPairs)),'Normalization','cdf_cleaned','DisplayStyle','stairs','EdgeColor',[1,0.5,0],'LineWidth',2);
hold on; histogram(histcounts( block1_cleaned.diffPairs,unique( block1_cleaned.diffPairs)),'Normalization','cdf_cleaned','FaceColor',[1,0.5,0],'LineWidth',2,'EdgeColor','none');
axis square;
xlabel('number of synaptic pairs');
ylabel('probability');
legend({'b1->b2','b1->b1'});

subplot(4,4,2)
histogram(histcounts(block2_cleaned.b2tob1diffPairs,unique( block2_cleaned.b2tob1diffPairs)),'Normalization','cdf_cleaned','DisplayStyle','stairs','EdgeColor',[0,0.5,1],'LineWidth',2);
hold on; histogram(histcounts( block2_cleaned.diffPairs,unique( block2_cleaned.diffPairs)),'Normalization','cdf_cleaned','FaceColor',[0,0.5,1],'LineWidth',2,'EdgeColor','none');
axis square;
legend({'b2->b1','b2->b2'});
%%

figure

% subplot(2,5,1)
% scatter(block1_cleaned.pairs_single(:,1),block1_cleaned.pairs_single(:,2),'MarkerFaceColor',[1,0.5,0],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
% hold on;
% f = showfit(ezfit(block1_cleaned.pairs_single(:,1),block1_cleaned.pairs_single(:,2),'affine'),'fitcolor','k','dispeqboxmode','off');
% %set(gca, 'XScale','log','YScale','log');
% axis square;
% text(10e2,10e3,sprintf('r = %0.2f',f.r));
% 

subplot(2,5,2)
scatter(block1_cleaned.dualpairs(:,1),block1_cleaned.dualpairs(:,2),'MarkerFaceColor',[1,0.5,0],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(block1_cleaned.dualpairs(:,1),block1_cleaned.dualpairs(:,2),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(10e2,10e3,sprintf('r = %0.2f',f.r));

subplot(2,5,3)
scatter(block1_cleaned.dualtriple(:,1),block1_cleaned.dualtriple(:,2),'MarkerFaceColor',[1,0.5,0],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(block1_cleaned.dualtriple(:,1),block1_cleaned.dualtriple(:,2),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(10e2,10e3,sprintf('r = %0.2f',f.r));

subplot(2,5,4)
scatter(block1_cleaned.dualquart(:,1),block1_cleaned.dualquart(:,2),'MarkerFaceColor',[1,0.5,0],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(block1_cleaned.dualquart(:,1),block1_cleaned.dualquart(:,2),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(10e2,10e3,sprintf('r = %0.2f',f.r));

subplot(2,5,5)
scatter(block1_cleaned.dualpent(:,1),block1_cleaned.dualpent(:,2),'MarkerFaceColor',[1,0.5,0],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(block1_cleaned.dualpent(:,1),block1_cleaned.dualpent(:,2),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(10e2,10e3,sprintf('r = %0.2f',f.r));



% subplot(2,5,6)
% scatter(block2_cleaned_gamma038.dualpairs(:,1),block2_cleaned_gamma038.dualpairs(:,2),'MarkerFaceColor',[0,0.5,1],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
% hold on;
% f = showfit(ezfit(block2_cleaned_gamma038.dualpairs(:,1),block2_cleaned_gamma038.dualpairs(:,2),'affine'),'fitcolor','k','dispeqboxmode','off');
% %set(gca, 'XScale','log','YScale','log');
% axis square;
% text(10e2,10e3,sprintf('r = %0.2f',f.r));


subplot(2,5,7)
scatter(block2_cleaned.dualpairs(:,1),block2_cleaned.dualpairs(:,2),'MarkerFaceColor',[0,0.5,1],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(block2_cleaned.dualpairs(:,1),block2_cleaned.dualpairs(:,2),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(10e2,10e3,sprintf('r = %0.2f',f.r));

subplot(2,5,8)
scatter(block2_cleaned.dualtriple(:,1),block2_cleaned.dualtriple(:,2),'MarkerFaceColor',[0,0.5,1],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(block2_cleaned.dualtriple(:,1),block2_cleaned.dualtriple(:,2),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(10e2,10e3,sprintf('r = %0.2f',f.r));

subplot(2,5,9)
scatter(block2_cleaned.dualquart(:,1),block2_cleaned.dualquart(:,2),'MarkerFaceColor',[0,0.5,1],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(block2_cleaned.dualquart(:,1),block2_cleaned.dualquart(:,2),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(10e2,10e3,sprintf('r = %0.2f',f.r));

subplot(2,5,10)
scatter(block2_cleaned.dualpent(:,1),block2_cleaned.dualpent(:,2),'MarkerFaceColor',[0,0.5,1],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(block2_cleaned.dualpent(:,1),block2_cleaned.dualpent(:,2),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(10e2,10e3,sprintf('r = %0.2f',f.r));


figure;


subplot(2,5,2)
scatter(log10(block1_cleaned.dualpairs(:,1)),log10(block1_cleaned.dualpairs(:,2)),'MarkerFaceColor',[1,0.5,0],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(log10(block1_cleaned.dualpairs(:,1)),log10(block1_cleaned.dualpairs(:,2)),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gc);
axis square;
text(4,4,sprintf('r = %0.2f',f.r));


subplot(2,5,3)
scatter(log10(block1_cleaned.dualtriple(:,1)),log10(block1_cleaned.dualtriple(:,2)),'MarkerFaceColor',[1,0.5,0],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(log10(block1_cleaned.dualtriple(:,1)),log10(block1_cleaned.dualtriple(:,2)),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(4,4,sprintf('r = %0.2f',f.r));

subplot(2,5,4)
scatter(log10(block1_cleaned.dualquart(:,1)),log10(block1_cleaned.dualquart(:,2)),'MarkerFaceColor',[1,0.5,0],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(log10(block1_cleaned.dualquart(:,1)),log10(block1_cleaned.dualquart(:,2)),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(4,4,sprintf('r = %0.2f',f.r));

subplot(2,5,5)
scatter(log10(block1_cleaned.dualpent(:,1)),log10(block1_cleaned.dualpent(:,2)),'MarkerFaceColor',[1,0.5,0],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(log10(block1_cleaned.dualpent(:,1)),log10(block1_cleaned.dualpent(:,2)),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(4,4,sprintf('r = %0.2f',f.r));



% subplot(2,5,6)
% scatter(block2_cleaned_gamma038.dualpairs(:,1),block2_cleaned_gamma038.dualpairs(:,2),'MarkerFaceColor',[0,0.5,1],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
% hold on;
% f = showfit(ezfit(block2_cleaned_gamma038.dualpairs(:,1),block2_cleaned_gamma038.dualpairs(:,2),'affine'),'fitcolor','k','dispeqboxmode','off');
% %set(gca, 'XScale','log','YScale','log');
% axis square;
% text(10e2,10e3,sprintf('r = %0.2f',f.r));


subplot(2,5,7)
scatter(log10(block2_cleaned.dualpairs(:,1)),log10(block2_cleaned.dualpairs(:,2)),'MarkerFaceColor',[0,0.5,1],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(log10(block2_cleaned.dualpairs(:,1)),log10(block2_cleaned.dualpairs(:,2)),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(4,4,sprintf('r = %0.2f',f.r));

subplot(2,5,8)
scatter(log10(block2_cleaned.dualtriple(:,1)),log10(block2_cleaned.dualtriple(:,2)),'MarkerFaceColor',[0,0.5,1],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(log10(block2_cleaned.dualtriple(:,1)),log10(block2_cleaned.dualtriple(:,2)),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(4,4,sprintf('r = %0.2f',f.r));

subplot(2,5,9)
scatter(log10(block2_cleaned.dualquart(:,1)),log10(block2_cleaned.dualquart(:,2)),'MarkerFaceColor',[0,0.5,1],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(log10(block2_cleaned.dualquart(:,1)),log10(block2_cleaned.dualquart(:,2)),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(4,4,sprintf('r = %0.2f',f.r));

subplot(2,5,10)
scatter(log10(block2_cleaned.dualpent(:,1)),log10(block2_cleaned.dualpent(:,2)),'MarkerFaceColor',[0,0.5,1],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(log10(block2_cleaned.dualpent(:,1)),log10(block2_cleaned.dualpent(:,2)),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(4,4,sprintf('r = %0.2f',f.r));

%%
figure;
subplot(4,4,1)
histogram(block1_cleaned.b1tob1size,'binWidth',500,'FaceColor',[1,0.5,0],'EdgeColor','none');
hold on
histogram(block2_cleaned.b2tob2size,'binWidth',500,'FaceColor',[0,0.5,1],'EdgeColor','none');
legend({'b1-->b1','b2-->b2'});
xlabel('voxels');
ylabel('count');
axis square;
box off;
offsetAxes(gca);

subplot(4,4,2)
histogram(log10(block1_cleaned.b1tob1size),50,'FaceColor',[1,0.5,0],'EdgeColor','none');
hold on
histogram(log10(block2_cleaned.b2tob2size),50,'FaceColor',[0,0.5,1],'EdgeColor','none');
legend({'b1-->b1','b2-->b2'});
xlabel('log voxels');
ylabel('count');
axis square;
box off;
set(gca,'XTickLabels',[10^2,10^3,10^4]);
offsetAxes(gca);

subplot(4,4,3)
histogram(block1_cleaned.b1tob2size,'binWidth',500,'EdgeColor',[1,0.5,0],'DisplayStyle','stairs','LineWidth',2);
hold on;
histogram(block2_cleaned.b2tob1size,'binWidth',500,'EdgeColor',[0,0.5,1],'DisplayStyle','stairs','LineWidth',2);
legend({'b1-->b2','b2-->b1'});
axis square;
box off;
offsetAxes(gca);

subplot(4,4,4)
histogram(log10(block1_cleaned.b1tob2size),20,'EdgeColor',[1,0.5,0],'DisplayStyle','stairs','LineWidth',2);
hold on;
histogram(log10(block2_cleaned.b2tob1size),20,'EdgeColor',[0,0.5,1],'DisplayStyle','stairs','LineWidth',2);
legend({'b1-->b2','b2-->b1'});
axis square;
xlabel('log voxels');
set(gca,'XTickLabels',[10^2,10^3,10^4]);
box off;
offsetAxes(gca);



subplot(4,4,5)
heatmap([mean(block1_cleaned.b1tob1size),mean(block2_cleaned.b2tob1size);mean(block1_cleaned.b1tob2size),mean(block2_cleaned.b2tob2size)],...
    'ColorbarVisible','off','XData',{'b1','b2'},'YData',{'b1','b2'},'FontSize',10);
title('mean synapse size (vx)');

%% size vs numbers

block1_cleaned.img_2 = PositionDensity(block1_cleaned.cellIDs,2,0);
block1_cleaned.img_3 = PositionDensity(block1_cleaned.cellIDs,3,0);
block1_cleaned.img_4 = PositionDensity(block1_cleaned.cellIDs,4,0);
block1_cleaned.img_5 = PositionDensity(block1_cleaned.cellIDs,5,0);

block2_cleaned.noVest = block2_cleaned.cellIDs(~isVestibular(block2_cleaned.cellIDs,df_cleaned));

block2_cleaned.img_2 = PositionDensity(block2_cleaned.noVest,2,0);
block2_cleaned.img_3 = PositionDensity(block2_cleaned.noVest,3,0);
block2_cleaned.img_4 = PositionDensity(block2_cleaned.noVest,4,0);
block2_cleaned.img_5 = PositionDensity(block2_cleaned.noVest,5,0);

% reshaped matrix

block1_cleaned.img_2reshape = reshape(cell2mat(block1_cleaned.img_2),[],1);
block2_cleaned.img_2reshape = reshape(cell2mat(block2_cleaned.img_2),[],1);

block1_cleaned.img_3reshape = reshape(cell2mat(block1_cleaned.img_3),[],1);
block2_cleaned.img_3reshape = reshape(cell2mat(block2_cleaned.img_3),[],1);

block1_cleaned.img_4reshape = reshape(cell2mat(block1_cleaned.img_4),[],1);
block2_cleaned.img_4reshape = reshape(cell2mat(block2_cleaned.img_4),[],1);

block1_cleaned.img_5reshape = reshape(cell2mat(block1_cleaned.img_5),[],1);
block2_cleaned.img_5reshape = reshape(cell2mat(block2_cleaned.img_5),[],1);

%
block1_cleaned.img_2_90 = sum(block1_cleaned.img_2reshape == 90);
block2_cleaned.img_2_90 = sum(block2_cleaned.img_2reshape == 90);

block1_cleaned.img_3_90 = sum(block1_cleaned.img_3reshape == 90);
block2_cleaned.img_3_90 = sum(block2_cleaned.img_3reshape == 90);

block1_cleaned.img_4_90 = sum(block1_cleaned.img_4reshape == 90);
block2_cleaned.img_4_90 = sum(block2_cleaned.img_4reshape == 90);

block1_cleaned.img_5_90 = sum(block1_cleaned.img_5reshape == 90);
block2_cleaned.img_5_90 = sum(block2_cleaned.img_5reshape == 90);
%%
% subplot(4,4,1)
% 
% plot([1,2,3,4],[block1_cleaned.img_2_90/numel(block1_cleaned.img_2reshape);...
%                 block1_cleaned.img_3_90/numel(block1_cleaned.img_3reshape);...
%                 block1_cleaned.img_4_90/numel(block1_cleaned.img_4reshape);...
%                 block1_cleaned.img_5_90/numel(block1_cleaned.img_5reshape)],'-ko');
% hold on;
% ylabel('Eyepos pixel fraction (mod1)');
% yyaxis right
% plot([1,2,3,4],[block2_cleaned_gamma038.img_2_90/numel(block2_cleaned_gamma038.img_2reshape);...
%                 block2_cleaned_gamma038.img_3_90/numel(block2_cleaned_gamma038.img_3reshape);...
%                 block2_cleaned_gamma038.img_4_90/numel(block2_cleaned_gamma038.img_4reshape);...
%                 block2_cleaned_gamma038.img_5_90/numel(block2_cleaned_gamma038.img_5reshape)],'-ro');
% set(gca,'YColor','r');
% ylabel('Eyepos pixel fraction (mod2)');
% set(gca,'XTickLabels',{'~6x6x8','~9x9x12','~12x12x16','~16x16x20'},'XTickLabelRotation',45);
% axis square;
% box off;

subplot(4,4,1)
plot([1,2,3,4],[block1_cleaned.img_2_90, block1_cleaned.img_3_90,...
    block1_cleaned.img_4_90,block1_cleaned.img_5_90],'-o','color',[1,0.5,0],'LineWidth',2);
hold on
plot([1,2,3,4],[block2_cleaned.img_2_90, block2_cleaned.img_3_90,...
    block2_cleaned.img_4_90,block2_cleaned.img_5_90],'-o','color',[0,0.5,1],'LineWidth',2);
set(gca,'XTickLabels',{'~6x6x8','~9x9x12','~12x12x16','~16x16x20'},'XTickLabelRotation',45);
legend({'mod1','mod2'});
ylabel('eyepos neurons');
 axis square;
 box off;



%% positionNess

optimalScale = 2;

%  cell1 = block1_cleaned.cellIDs;
%  cell2 = block2_cleaned.noVest;

 cell1 = setdiff(block1_cleaned.cellIDs,block1_cleaned.noSoma);
 cell2 = setdiff(block2_cleaned.cellIDs,block2_cleaned.noSoma);
%control
% block1_cleaned.img_2 = PositionDensity(setdiff(block1_cleaned.cellIDs,block1_cleaned.noSoma),optimalScale,0);
% block2_cleaned.img_2 = PositionDensity(setdiff(block2_cleaned.cellIDs,block2_cleaned.noSoma),optimalScale,0);

block1_cleaned.img_2 = PositionDensity(cell1,optimalScale,0);
block2_cleaned.img_2 = PositionDensity(cell2,optimalScale,0);


block1_cleaned.img_2reshape = reshape(cell2mat(block1_cleaned.img_2),[],1);
block2_cleaned.img_2reshape = reshape(cell2mat(block2_cleaned.img_2),[],1);

nreps = 10;

for i = 1:nreps

% jitter 5um
block1_cleaned.img_2_5 = PositionDensity(cell1,optimalScale,5);
block2_cleaned.img_2_5 = PositionDensity(cell2,optimalScale,5);

block1_cleaned.img_2_5reshape = reshape(cell2mat(block1_cleaned.img_2_5),[],1);
block2_cleaned.img_2_5reshape = reshape(cell2mat(block2_cleaned.img_2_5),[],1);

block1_cleaned.posNeurons_5(i) = sum(block1_cleaned.img_2_5reshape  == 90);
block2_cleaned.posNeurons_5(i) = sum(block2_cleaned.img_2_5reshape  == 90);

% jitter 10um

block1_cleaned.img_2_10 = PositionDensity(cell1,optimalScale,10);
block2_cleaned.img_2_10 = PositionDensity(cell2,optimalScale,10);

block1_cleaned.img_2_10reshape = reshape(cell2mat(block1_cleaned.img_2_10),[],1);
block2_cleaned.img_2_10reshape = reshape(cell2mat(block2_cleaned.img_2_10),[],1);

block1_cleaned.posNeurons_10(i) = sum(block1_cleaned.img_2_10reshape  == 90);
block2_cleaned.posNeurons_10(i) = sum(block2_cleaned.img_2_10reshape  == 90);


% jitter 5um

block1_cleaned.img_2_15 = PositionDensity(cell1,optimalScale,15);
block2_cleaned.img_2_15 = PositionDensity(cell2,optimalScale,15);

block1_cleaned.img_2_15reshape = reshape(cell2mat(block1_cleaned.img_2_15),[],1);
block2_cleaned.img_2_15reshape = reshape(cell2mat(block2_cleaned.img_2_15),[],1);

block1_cleaned.posNeurons_15(i) = sum(block1_cleaned.img_2_15reshape  == 90);
block2_cleaned.posNeurons_15(i) = sum(block2_cleaned.img_2_15reshape  == 90);

end

%%
figure;
subplot(442)

errorbar([sum(block1_cleaned.img_2reshape == 90);...
           mean(block1_cleaned.posNeurons_5);...
            mean(block1_cleaned.posNeurons_10);...
             mean(block1_cleaned.posNeurons_15)],... 
         [0;std(block1_cleaned.posNeurons_5)/sqrt(20);...
             std(block1_cleaned.posNeurons_5)/sqrt(20);...
             std(block1_cleaned.posNeurons_5)/sqrt(20)],'-o','color',[1,0.5,0],'LineWidth',2);
hold on
errorbar([sum(block2_cleaned.img_2reshape == 90);...
           mean(block2_cleaned.posNeurons_5);...
            mean(block2_cleaned.posNeurons_10);...
             mean(block2_cleaned.posNeurons_15)],... 
         [0;std(block2_cleaned.posNeurons_5)/sqrt(20);...
             std(block2_cleaned.posNeurons_5)/sqrt(20);...
             std(block2_cleaned.posNeurons_5)/sqrt(20)],'-o','color',[0,0.5,1],'LineWidth',2);
         
         set(gca,'XTickLabels',{'centered','5\mu','10\mu','15\mu'},'XTickLabelRotation',35);
box off;
title('@ ~6x6x8');
axis square;
%offsetAxes(gca);


subplot(4,4,3)        
hold on;
errorbar([sum(block2_cleaned.img_2reshape == 90)/sum(block1_cleaned.img_2reshape == 90);...
                mean((block2_cleaned.posNeurons_5./block1_cleaned.posNeurons_5));...
                mean((block2_cleaned.posNeurons_10./block1_cleaned.posNeurons_10));...
                mean((block2_cleaned.posNeurons_15./block1_cleaned.posNeurons_15))],...
                [0;...
                std((block2_cleaned.posNeurons_5./block1_cleaned.posNeurons_5))/sqrt(20);...
                std((block2_cleaned.posNeurons_10./block1_cleaned.posNeurons_10))/sqrt(20);...
                std((block2_cleaned.posNeurons_15./block1_cleaned.posNeurons_15))/sqrt(20)],'-or');
set(gca,'XTickLabels',{'centered','5\mu','10\mu','15\mu'},'XTickLabelRotation',35);
ylabel('mod2/mod1');
box off;
title('@ ~6x6x8');
axis square;
%offsetAxes(gca);

%%
figure;
subplot(4,4,2)
plot([1,2,3,4],[sum(block1_cleaned.img_2reshape == 90);...
                sum(block1_cleaned.img_2_5reshape  == 90);...
                sum(block1_cleaned.img_2_10reshape ==90);...
                sum(block1_cleaned.img_2_15reshape ==90)],'-o','color',[1,0.5,0],'LineWidth',2);    
ylabel('eyepos neurons')
hold on;
plot([1,2,3,4],[sum(block2_cleaned.img_2reshape == 90);...
                sum(block2_cleaned.img_2_5reshape  == 90);...
                sum(block2_cleaned.img_2_10reshape ==90);...
                sum(block2_cleaned.img_2_15reshape ==90)],'-o','color',[0,0.5,1],'LineWidth',2);
set(gca,'XTickLabels',{'centered','5\mu','10\mu','15\mu'},'XTickLabelRotation',35);
box off;
title('@ ~6x6x8');
axis square;

subplot(4,4,3)        
hold on;
plot([1,2,3,4],[sum(block2_cleaned.img_2reshape == 90)/sum(block1_cleaned.img_2reshape == 90);...
                sum(block2_cleaned.img_2_5reshape  == 90)/sum(block1_cleaned.img_2_5reshape  == 90);...
                sum(block2_cleaned.img_2_10reshape ==90)/sum(block1_cleaned.img_2_10reshape ==90);...
                sum(block2_cleaned.img_2_15reshape ==90)/ sum(block1_cleaned.img_2_15reshape ==90)],'-or');
set(gca,'XTickLabels',{'centered','5\mu','10\mu','15\mu'},'XTickLabelRotation',35);
ylabel('mod2/mod1');
box off;
title('@ ~6x6x8');

axis square;

%%
subplot(4,4,2)

plot([1,2,3,4],[sum(block1_cleaned.img_2reshape == 90)/numel(block1_cleaned.img_2reshape);...
                sum(block1_cleaned.img_2_5reshape  == 90)/numel(block1_cleaned.img_2_5reshape);...
                sum(block1_cleaned.img_2_10reshape ==90)/ numel(block1_cleaned.img_2_10reshape);...
                sum(block1_cleaned.img_2_15reshape ==90)/numel(block1_cleaned.img_2_15reshape)],'-ok');          
ylabel('eyepos pixel fraction');
yyaxis right;

plot([1,2,3,4],[sum(block2_cleaned.img_2reshape == 90)/numel(block2_cleaned.img_2reshape);...
                sum(block2_cleaned.img_2_5reshape  == 90)/numel(block2_cleaned.img_2_5reshape);...
                sum(block2_cleaned.img_2_10reshape ==90)/ numel(block2_cleaned.img_2_10reshape);...
                sum(block2_cleaned.img_2_15reshape ==90)/numel(block2_cleaned.img_2_15reshape)],'-or');
set(gca,'YColor','r');
set(gca,'XTickLabels',{'centered','5\mu','10\mu','15\mu'},'XTickLabelRotation',35);
box off;
title('@ ~6x6x8');

subplot(4,4,3)

plot([1,2,3,4],[(sum(block2_cleaned.img_2reshape == 90)/numel(block2_cleaned.img_2reshape))/...
                (sum(block1_cleaned.img_2reshape == 90)/numel(block1_cleaned.img_2reshape));...
                 (sum(block2_cleaned.img_2_5reshape  == 90)/numel(block2_cleaned.img_2_5reshape))/...
                 (sum(block1_cleaned.img_2_5reshape  == 90)/numel(block1_cleaned.img_2_5reshape));...
                 (sum(block2_cleaned.img_2_10reshape ==90)/ numel(block2_cleaned.img_2_10reshape))/...
                 (sum(block1_cleaned.img_2_10reshape ==90)/ numel(block1_cleaned.img_2_10reshape));...
                 (sum(block2_cleaned.img_2_15reshape ==90)/numel(block2_cleaned.img_2_15reshape))/...
                 (sum(block1_cleaned.img_2_15reshape ==90)/numel(block1_cleaned.img_2_15reshape))],'-o');
set(gca,'XTickLabels',{'centered','5\mu','10\mu','15\mu'},'XTickLabelRotation',35);
box off;
ylabel('mod2/mod1');

%%
% 
% patchSize = [2*2*2*0.798,2*2*2*0.798,2/2*2*2;...
%              2*3*2*0.798,2*3*2*0.798,3/2*2*2;...
%              2*4*2*0.798,2*4*2*0.798,4/2*2*2;...
%              2*5*2*0.798,2*5*2*0.798,5/2*2*2];

% subplot(4,4,1)
% plot([1,2,3,4],[length(block2_cleaned_gamma038.img_2(block2_cleaned_gamma038.img_2 == 90))/length(block1_cleaned.img_2(block1_cleaned.img_2 == 90)),...
%    length(block2_cleaned_gamma038.img_3(block2_cleaned_gamma038.img_3 == 90))/length(block1_cleaned.img_3(block1_cleaned.img_3 == 90)),...
%    length(block2_cleaned_gamma038.img_4(block2_cleaned_gamma038.img_4 == 90))/length(block1_cleaned.img_4(block1_cleaned.img_4 == 90)),...
%    length(block2_cleaned_gamma038.img_5(block2_cleaned_gamma038.img_5 == 90))/length(block1_cleaned.img_5(block1_cleaned.img_5 == 90))],'-o');
% set(gca,'XTickLabels',{'~6x6x8','~9x9x12','~12x12x16','~16x16x20'},'XTickLabelRotation',45);
% box off;
% daspect([2,2,1]);
% xlabel('patch size (\mum)')
% ylabel('mod2/mod1 cells')

subplot(4,4,2)
plot([1,2,3,4],...
    [length(block2_cleaned.img_2(block2_cleaned.img_2 == 90))/length(block1_cleaned.img_2(block1_cleaned.img_2 == 90)),...
    length(block2_cleaned.img_2_5(block2_cleaned.img_2_5 == 90))/length(block1_cleaned.img_2_5(block1_cleaned.img_2_5 == 90)),...
    length(block2_cleaned.img_2_10(block2_cleaned.img_2_10 == 90))/length(block1_cleaned.img_2_10(block1_cleaned.img_2_10 == 90)),...
    length(block2_cleaned.img_2_15(block2_cleaned.img_2_15 == 90))/length(block1_cleaned.img_2(block1_cleaned.img_2_15 == 90))],...
    '-o');
set(gca,'XTickLabels',{'centered','5\mu','10\mu','15\mu'},'XTickLabelRotation',35);
box off;
daspect([1,2,1]);
xlabel('centroid location')
ylabel('mod2/mod1 cells')
title('@ ~6x6x8')

subplot(4,4,3)

plot([1,2,3,4],...
    [length(block1_cleaned.img_2(block1_cleaned.img_2 == 90)),...
    length(block1_cleaned.img_2_5(block1_cleaned.img_2_5 == 90)),...
    length(block1_cleaned.img_2_10(block1_cleaned.img_2_10 == 90)),...
    length(block1_cleaned.img_2(block1_cleaned.img_2_15 == 90))],...
    '-o');
hold on;
plot([1,2,3,4],...
    [length(block2_cleaned.img_2(block2_cleaned.img_2 == 90)),...
    length(block2_cleaned.img_2_5(block2_cleaned.img_2_5 == 90)),...
    length(block2_cleaned.img_2_10(block2_cleaned.img_2_10 == 90)),...
    length(block2_cleaned.img_2_15(block2_cleaned.img_2_15 == 90))],...
        '-o');    
set(gca,'XTickLabels',{'centered','5\mu','10\mu','15\mu'},'XTickLabelRotation',35);
box off;
daspect([1,4,1]);
xlabel('centroid location')
ylabel('number of cells')
legend({'mod1','mod2'})


% block1_cleaned_gamma09_WithNull.img = PositionDensity(block1_cleaned_gamma09_WithNull.cellIDs);
% block2_cleaned_gamma038_gamma09_WithNull.img = PositionDensity(block2_cleaned_gamma038_gamma09_WithNull.cellIDs);
% block3_gamma09_WithNull.img = PositionDensity(block3_gamma09_WithNull.cellIDs);
% 


figure(1)
for i = 1:length(block1_cleaned.cellIDs)
    subplot(13,13,i);
    imagesc(block1_cleaned.img_2(:,:,i));
    axis off;
    axis square;
    box on;
    title(block1_cleaned.cellIDs(i),'FontSize',8);
end
sgtitle('block1_cleaned');
%colormap(colorcet('L19'));

figure(2);
for i = 1:length(block2_cleaned.cellIDs)
    subplot(13,13,i);
    imagesc(block2_cleaned.img_2(:,:,i));
    axis off;
    axis square;
    box on;
    title(block2_cleaned.cellIDs(i),'FontSize',8);
end
sgtitle('block2_cleaned_gamma038');
% 
% figure(3);
% for i = 1:length(block1_cleaned_gamma09_WithNull.cellIDs)
%     subplot(13,13,i);
%     imagesc(block1_cleaned_gamma09_WithNull.img(:,:,i));
%     axis off;
%     axis square;
%     box on;
%     title(block1_cleaned_gamma09_WithNull.cellIDs(i),'FontSize',8);
% end
% sgtitle('block1_cleaned_gamma09_WithNull');
% 
% figure(4);
% for i = 1:length(block2_cleaned_gamma038_gamma09_WithNull.cellIDs)
%     subplot(13,13,i);
%     imagesc(block2_cleaned_gamma038_gamma09_WithNull.img(:,:,i));
%     axis off;
%     axis square;
%     box on;
%     title(block2_cleaned_gamma038_gamma09_WithNull.cellIDs(i),'FontSize',8);
% end
% sgtitle('block2_cleaned_gamma038_gamma09_WithNull');
% 
% figure(5);
% for i = 1:length(block3_gamma09_WithNull.cellIDs)
%     subplot(13,13,i);
%     imagesc(block3_gamma09_WithNull.img(:,:,i));
%     axis off;
%     axis square;
%     box on;
%     title(block3_gamma09_WithNull.cellIDs(i),'FontSize',8);
% end
% sgtitle('block3_gamma09_WithNull');

%%

figure;

subplot(1,2,1)
polarhistogram(block1_cleaned.img_2(:),25,'Facecolor',[1,0.5,0])
hold on
polarhistogram(block2_cleaned.img_2(:),25,'Facecolor',[0,0.5,1])


subplot(1,2,2)
polarhistogram(block1_cleaned.img_2_10(:),25,'Facecolor',[1,0.5,0])
hold on
polarhistogram(block2_cleaned.img_2_10(:),25,'Facecolor',[0,0.5,1])



%%

figure;

subplot(4,4,2)
histogram(block1_cleaned.img_2,'BinWidth',5,'EdgeColor',[1,0.5,0],'DisplayStyle','stairs','LineWidth',2);
hold on
histogram(block2_cleaned.img_2,'BinWidth',5,'EdgeColor',[0,0.5,1],'DisplayStyle','stairs','LineWidth',2);
axis square;
box off;
xlabel('grey value');
ylabel('freq');


subplot(4,4,5)
histogram(block1_cleaned.img_2_5,'BinWidth',5,'EdgeColor',[1,0.5,0],'DisplayStyle','stairs','LineWidth',2);
hold on
histogram(block2_cleaned.img_2_5,'BinWidth',5,'EdgeColor',[0,0.5,1],'DisplayStyle','stairs','LineWidth',2);
%set(gca,'YLim',[0,0.04]);
axis square;
box off; 
title('5\mum jitter');

subplot(4,4,6)
histogram(block1_cleaned.img_2_10,'BinWidth',5,'EdgeColor',[1,0.5,0],'DisplayStyle','stairs','LineWidth',2);
hold on
histogram(block2_cleaned.img_2_10,'BinWidth',5,'EdgeColor',[0,0.5,1],'DisplayStyle','stairs','LineWidth',2);
%set(gca,'YLim',[0,0.04]);
axis square;
box off;
title('10\mum jitter');


subplot(4,4,7)
histogram(block1_cleaned.img_2_15,'BinWidth',5,'EdgeColor',[1,0.5,0],'DisplayStyle','stairs','LineWidth',2);
hold on
histogram(block2_cleaned.img_2_15,'BinWidth',5,'EdgeColor',[0,0.5,1],'DisplayStyle','stairs','LineWidth',2);
%set(gca,'YLim',[0,0.04]);
axis square;
box off;
title('15\mum jitter');



%%

errorbar([nanmean(block1_cleaned.img{2}(:));nanmean(block1_cleaned.img{4}(:));nanmean(block1_cleaned.img{6}(:));nanmean(block1_cleaned.img{8}(:));nanmean(block1_cleaned.img{10}(:));nanmean(block1_cleaned.img{12}(:))],...
    [nanstd(block1_cleaned.img{2}(:));nanstd(block1_cleaned.img{4}(:));nanstd(block1_cleaned.img{6}(:));nanstd(block1_cleaned.img{8}(:));nanstd(block1_cleaned.img{10}(:));nanstd(block1_cleaned.img{12}(:))]./[157;157;157;157;157;157])

hold on;
errorbar([nanmean(block2_cleaned.img{2}(:));nanmean(block2_cleaned.img{4}(:));nanmean(block2_cleaned.img{6}(:));nanmean(block2_cleaned.img{8}(:));nanmean(block2_cleaned.img{10}(:));nanmean(block2_cleaned.img{12}(:))],...
    [nanstd(block2_cleaned.img{2}(:));nanstd(block2_cleaned.img{4}(:));nanstd(block2_cleaned.img{6}(:));nanstd(block2_cleaned.img{8}(:));nanstd(block2_cleaned.img{10}(:));nanstd(block2_cleaned.img{12}(:))]./[116;116;116;116;116;116])
