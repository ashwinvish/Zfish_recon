clc
clear 
%
load AllCells.mat
load df_cleaned.mat % need access since large file
load ConnMatrixPre_cleaned.mat

allAxonID = df_cleaned.presyn_segid(:);
allAxonID_recon =  allAxonID(allAxonID<1e5);
allAxonSites = [df_cleaned.presyn_x(:),df_cleaned.presyn_y(:),df_cleaned.presyn_z(:)];
allAxonSites_recon = allAxonSites(allAxonID<1e5,:);

[a,b] = ismember(allAxonID_recon,df_cleaned.presyn_segid(:));
allPSDsize_recon = df_cleaned.size(a);

allAxonSites_recon = TransformPoints(allAxonSites_recon,0); % covert from voxel space to micron space


actual_synapses = zeros(size(ConnMatrixPre_cleaned));
potentialSynapses_2 = zeros(size(ConnMatrixPre_cleaned));
potentialSynapses_5 = zeros(size(ConnMatrixPre_cleaned));
potentialSynapses_10 = zeros(size(ConnMatrixPre_cleaned));

actual_synSize = [];
Pot2_SynSize = [];
Pot5_SynSize = [];
Pot10_SynSize = [];

for i = 1:length(AllCells)
    [PrePartners,prePSD] = SynapticPartners(AllCells(i),1,df_cleaned);
    inputs = PrePartnerCoordinates(prePSD,df_cleaned); % location of presynaptic partner, on the axon
    
    if ~isempty(inputs)
        inputs = TransformPoints(inputs,0);
        temp2 = zeros(size(inputs,1),length(AllCells));
        temp5 = zeros(size(inputs,1),length(AllCells));
        temp10 = zeros(size(inputs,1),length(AllCells));

        for j = 1:size(inputs,1)
            distMat = pdist2(allAxonSites_recon,inputs(j,:));
            
            actual_syn = allAxonID_recon(distMat == 0);
            PotSyn_2 = allAxonID_recon(distMat <= 2);
            PotSyn_5 = allAxonID_recon(distMat <= 5);
            PotSyn_10 = allAxonID_recon(distMat <= 10);
            
                        
            [n0,e0] = edgeCounts(actual_syn,AllCells);
            temp0(j,e0) = n0;


            [n2,e2] = edgeCounts(PotSyn_2,AllCells);
            temp2(j,e2) = n2;
                
            [n5,e5] = edgeCounts(PotSyn_5,AllCells);
            temp5(j,e5) = n5;
            
            [n10,e10] = edgeCounts(PotSyn_10,AllCells);
            temp10(j,e10) = n10;
            
            clear actual_syn;
            clear PotSyn_2;
            clear PotSyn_5;
            clear PotSyn_10;
            
        end
        
    potentialSynapses_2(i,:) = sum(temp2,1);
    potentialSynapses_5(i,:) = sum(temp5,1);
    potentialSynapses_10(i,:) = sum(temp10,1);
       
    
    clear PrePartners
    clear prePSD
    clear inputs
    clear distMat
    end
end

load block2_cleaned_gamma038_08072020.mat
load block1_cleaned_gamma038_08072020.mat

ABDr_CellIDs = [77305,77672,82194,77300,77709,77710,77705,82146,77301,77302,82143,77648,82145,77661,82140];
ABDc_CellIDs = [77628,77292,77296,77688,77154,77658,82195,82212,77295,77652,77646,77682,82197,81172,77657,77654,82213];
ABDIr_CellIDs = [78556,78552,77618,77668,77158,77665,77634,78547,77150 77631,77886,78553];
ABDIc_CellIDs = [79066,77640,77692,77625,77643,77144,77148,79051,77641 78574];

allABD = [ABDr_CellIDs,ABDc_CellIDs,ABDIr_CellIDs,ABDIc_CellIDs];

RScells_Large = [76202,77267,77268,76562,77931];
RScells_Small = [77260 76562 77259 77265 77267 77268 ...
    77694 77695 77931 78566 78577 78940 82217 82220];

[~,matIndex_1] = ismember(block1_cleaned.cellIDs,AllCells);
[~,matIndex_2] = ismember(block2_cleaned.cellIDs,AllCells);
[~,matIndex_ABD] = ismember(allABD,AllCells);
[~,matIndex_RSl] = ismember(RScells_Large,AllCells);


matIndex_2_ABD = [matIndex_2,matIndex_ABD];


PotentialConnectome_22_2 = potentialSynapses_2(matIndex_2_ABD,matIndex_2_ABD);
PotentialConnectome_22_2 = PotentialConnectome_22_2 - diag(diag(PotentialConnectome_22_2));
PotentialConnectome_22_2(:,size(block2_cleaned.cellIDs,2):end) = 0;

PotentialConnectome_22_5 = potentialSynapses_5(matIndex_2_ABD,matIndex_2_ABD);
PotentialConnectome_22_5 = PotentialConnectome_22_5 - diag(diag(PotentialConnectome_22_5));
PotentialConnectome_22_5(:,size(block2_cleaned.cellIDs,2):end) = 0;


PotentialConnectome_22_10 = potentialSynapses_10(matIndex_2_ABD,matIndex_2_ABD);
PotentialConnectome_22_10 = PotentialConnectome_22_10 - diag(diag(PotentialConnectome_22_10));
PotentialConnectome_22_10(:,size(block2_cleaned.cellIDs,2):end) = 0;



for i = 1:length([matIndex_2,matIndex_ABD])
    TotalSynapsesPot_2(i) = sum(potentialSynapses_2(matIndex_2_ABD(i),:));
    TotalSynapsesPot_5(i) = sum(potentialSynapses_5(matIndex_2_ABD(i),:));
    TotalSynapsesPot_10(i) = sum(potentialSynapses_10(matIndex_2_ABD(i),:));
end



% save('TotalSynapsesPot_2.mat','TotalSynapsesPot_2');
% save('TotalSynapsesPot_5.mat','TotalSynapsesPot_5');
% save('TotalSynapsesPot_10.mat','TotalSynapsesPot_10');
% 
% save('PotentialConnectome_22_2.mat','PotentialConnectome_22_2');
% save('PotentialConnectome_22_5.mat','PotentialConnectome_22_5');
% save('PotentialConnectome_22_10.mat','PotentialConnectome_22_10');
% 
% save('potentialSynapses_2.mat','potentialSynapses_2');
% save('potentialSynapses_5.mat','potentialSynapses_5');
% save('potentialSynapses_10.mat','potentialSynapses_10');


%% potential synapses from i to j

load block2_cleaned_gamma038_08072020.mat
load block1_cleaned_gamma038_08072020.mat

load potentialSynapses_2.mat
load potentialSynapses_5.mat
load potentialSynapses_10.mat

load AllCells.mat
load ConnMatrixPre_cleaned.mat
load df_cleaned.mat



ABDr_CellIDs = [77305,77672,82194,77300,77709,77710,77705,82146,77301,77302,82143,77648,82145,77661,82140];
ABDc_CellIDs = [77628,77292,77296,77688,77154,77658,82195,82212,77295,77652,77646,77682,82197,81172,77657,77654,82213];
ABDIr_CellIDs = [78556,78552,77618,77668,77158,77665,77634,78547,77150 77631,77886,78553];
ABDIc_CellIDs = [79066,77640,77692,77625,77643,77144,77148,79051,77641 78574];

allABD = [ABDr_CellIDs,ABDc_CellIDs,ABDIr_CellIDs,ABDIc_CellIDs];

RScells_Large = [76202,77267,77268,76562,77931];
RScells_Small = [77260 76562 77259 77265 77267 77268 ...
    77694 77695 77931 78566 78577 78940 82217 82220];

[~,matIndex_1] = ismember(block1_cleaned.cellIDs,AllCells);
[~,matIndex_2] = ismember(block2_cleaned.cellIDs,AllCells);
[~,matIndex_ABD] = ismember(allABD,AllCells);
[~,matIndex_RSl] = ismember(RScells_Large,AllCells);


matIndex_2_ABD = [matIndex_2,matIndex_ABD];



block1.TrueSynapse(1,1) = sum(sum(ConnMatrixPre_cleaned(matIndex_1,matIndex_1)));
block1.TrueSynapse(1,2) = sum(sum(ConnMatrixPre_cleaned(matIndex_2,matIndex_1)));
block1.TrueSynapse(2,2) = sum(sum(ConnMatrixPre_cleaned(matIndex_2,matIndex_2)));
block1.TrueSynapse(2,1) = sum(sum(ConnMatrixPre_cleaned(matIndex_1,matIndex_2)));
block1.Denom = [size(block1_cleaned.cellIDs,2)*size(block1_cleaned.cellIDs,2),...
    size(block1_cleaned.cellIDs,2)*size(block2_cleaned.cellIDs,2);...
    size(block1_cleaned.cellIDs,2)*size(block2_cleaned.cellIDs,2),...
    size(block2_cleaned.cellIDs,2)*size(block2_cleaned.cellIDs,2)];

block1.PotentialSynapse_2(1,1) = sum(sum(potentialSynapses_2(matIndex_1,matIndex_1)));
block1.PotentialSynapse_2(1,2) = sum(sum(potentialSynapses_2(matIndex_2,matIndex_1)));
block1.PotentialSynapse_2(2,2) = sum(sum(potentialSynapses_2(matIndex_2,matIndex_2)));
block1.PotentialSynapse_2(2,1) = sum(sum(potentialSynapses_2(matIndex_1,matIndex_2)));

block1.PotentialSynapse_5(1,1) = sum(sum(potentialSynapses_5(matIndex_1,matIndex_1)));
block1.PotentialSynapse_5(1,2) = sum(sum(potentialSynapses_5(matIndex_2,matIndex_1)));
block1.PotentialSynapse_5(2,2) = sum(sum(potentialSynapses_5(matIndex_2,matIndex_2)));
block1.PotentialSynapse_5(2,1) = sum(sum(potentialSynapses_5(matIndex_1,matIndex_2)));

block1.PotentialSynapse_10(1,1) = sum(sum(potentialSynapses_10(matIndex_1,matIndex_1)));
block1.PotentialSynapse_10(1,2) = sum(sum(potentialSynapses_10(matIndex_2,matIndex_1)));
block1.PotentialSynapse_10(2,2) = sum(sum(potentialSynapses_10(matIndex_2,matIndex_2)));
block1.PotentialSynapse_10(2,1) = sum(sum(potentialSynapses_10(matIndex_1,matIndex_2)));

% block1.b1tob1size = synSize(block1_cleaned.cellIDs,block1_cleaned.cellIDs,df_cleaned);
% block1.b1tob2size = synSize(block1_cleaned.cellIDs,block2_cleaned.cellIDs,df_cleaned);
% block1.b1toABDsize = synSize(block1_cleaned.cellIDs,allABD,df_cleaned);
% block1.b1toRSsize = synSize(block1_cleaned.cellIDs,RScells_Large,df_cleaned);
% 
% block2.b2tob1size = synSize(block2_cleaned.cellIDs,block1_cleaned.cellIDs,df_cleaned);
% block2.b2tob2size = synSize(block2_cleaned.cellIDs,block2_cleaned.cellIDs,df_cleaned);
% block2.b2toABDsize = synSize(block2_cleaned.cellIDs,allABD,df_cleaned);
% block2.b2toRSsize = synSize(block2_cleaned.cellIDs,RScells_Large,df_cleaned);
% 

figure;
subplot(4,4,1)
plot([1,2,3,4],[DiagRatio(block1.TrueSynapse./block1.Denom ),DiagRatio(block1.PotentialSynapse_2./block1.Denom ),...
    DiagRatio(block1.PotentialSynapse_5./block1.Denom ),DiagRatio(block1.PotentialSynapse_10./block1.Denom)],'-ko','LineWidth',2);
set(gca,'XTick',[1:4],'XTickLabels',[{'true'},2,5,10],'YLim',[0,8]);
box off;
daspect([1,1,1]);
%offsetAxes(gca);
xlabel('radius (\mum)');
ylabel('Sum within-module/Sum between-module')


block1.ABDSynapse(1,1) = sum(sum(ConnMatrixPre_cleaned(matIndex_ABD,matIndex_1)));
block1.ABDSynapse(1,2) = sum(sum(potentialSynapses_2(matIndex_ABD,matIndex_1)));
block1.ABDSynapse(1,3) = sum(sum(potentialSynapses_5(matIndex_ABD,matIndex_1)));
block1.ABDSynapse(1,4) = sum(sum(potentialSynapses_10(matIndex_ABD,matIndex_1)));
block1.ABDdenom = (size(block1_cleaned.cellIDs,2)*size(allABD,2));

block2.ABDSynapse(1,1) = sum(sum(ConnMatrixPre_cleaned(matIndex_ABD,matIndex_2)));
block2.ABDSynapse(1,2) = sum(sum(potentialSynapses_2(matIndex_ABD,matIndex_2)));
block2.ABDSynapse(1,3) = sum(sum(potentialSynapses_5(matIndex_ABD,matIndex_2)));
block2.ABDSynapse(1,4) = sum(sum(potentialSynapses_10(matIndex_ABD,matIndex_2)));
block2.ABDdenom = (size(block2_cleaned.cellIDs,2)*size(allABD,2));





block1.RSlSynapse(1,1) = sum(sum(ConnMatrixPre_cleaned(matIndex_RSl,matIndex_1)));
block1.RSlSynapse(1,2) = sum(sum(potentialSynapses_2(matIndex_RSl,matIndex_1)));
block1.RSlSynapse(1,3) = sum(sum(potentialSynapses_5(matIndex_RSl,matIndex_1)));
block1.RSlSynapse(1,4) = sum(sum(potentialSynapses_10(matIndex_RSl,matIndex_1)));
block1.RSdenom = (size(block1_cleaned.cellIDs,2)*size(RScells_Large,2));


block2.RSlSynapse(1,1) = sum(sum(ConnMatrixPre_cleaned(matIndex_RSl,matIndex_2)));
block2.RSlSynapse(1,2) = sum(sum(potentialSynapses_2(matIndex_RSl,matIndex_2)));
block2.RSlSynapse(1,3) = sum(sum(potentialSynapses_5(matIndex_RSl,matIndex_2)));
block2.RSlSynapse(1,4) = sum(sum(potentialSynapses_10(matIndex_RSl,matIndex_2)));
block2.RSdenom = (size(block2_cleaned.cellIDs,2)*size(RScells_Large,2));


center2periphery_Area_norm = [block1.ABDSynapse(1,1)./(size(block1_cleaned.cellIDs,2)*size(allABD,2)),...
                         block2.ABDSynapse(1,1)./(size(block2_cleaned.cellIDs,2)*size(allABD,2));...
                         block1.RSlSynapse(1,1)./(size(block1_cleaned.cellIDs,2)*size(RScells_Large,2)),...
                         block2.RSlSynapse(1,1)./(size(block2_cleaned.cellIDs,2)*size(RScells_Large,2))];
                     
 center2periphery_Area = [block1.ABDSynapse(1,1),...
     block2.ABDSynapse(1,1);...
     block1.RSlSynapse(1,1),...
     block2.RSlSynapse(1,1)];


subplot(4,4,2)
plot([1,2,3,4],(block2.ABDSynapse./block2.ABDdenom) ./ (block1.ABDSynapse./block1.ABDdenom),'-ko','LineWidth',2);
set(gca,'YLim',[0,18]);
xlabel('radius (\mum)');
ylabel('modO-->ABD/modA -->ABD')
box off;
yyaxis right
plot([1,2,3,4],(block1.RSlSynapse./block1.RSdenom)./ (block2.RSlSynapse./block2.RSdenom),'-ro','LineWidth',2);
set(gca,'XTick',[1:4],'XTickLabels',[{'true'},2,5,10]);
ylabel('modA-->RS/modO -->RS')
box off;
daspect([1,2,1]);
%offsetAxes(gca);

subplot(4,4,3)
plot([1,2,3,4],(block2.ABDSynapse./block2.ABDdenom) ./ (block2.RSlSynapse./block2.RSdenom),'-ko','LineWidth',2);
xlabel('radius (\mum)');
ylabel('modO-->ABD/modO -->RS')
box off;
yyaxis right
plot([1,2,3,4],(block1.RSlSynapse./block1.RSdenom)./ (block1.ABDSynapse./block1.ABDdenom),'-ro','LineWidth',2);
set(gca,'XTick',[1:4],'XTickLabels',[{'true'},2,5,10]);
ylabel('modA-->RS/modA -->ABD')
box off;


subplot(4,4,4)
%hold on;
plot([1,2,3,4],((block2.ABDSynapse./block2.ABDdenom) + (block1.RSlSynapse./block1.RSdenom)) ./ ((block2.RSlSynapse./ block2.RSdenom) + (block1.ABDSynapse./block1.ABDdenom)),...
    '-ko','LineWidth',2);
set(gca,'XTick',[1:4],'XTickLabels',[{'true'},2,5,10],'YLim',[0,12],'YTick',0:2:12,'YTickLabels',[0,2,4,6,8,10,12]);
ylabel('$\displaystyle\frac{Sum to preferred partner}{ sum to non-preferred partner}$','Interpreter','latex')
box off;
daspect([1,2,1])
%offsetAxes(gca);

%% submods


load block1_cleaned_submod_gamma038_08072020.mat
load block2_cleaned_submod_gamma038_08072020.mat

[~,matIndex_sub1] = ismember(block1_cleaned_submod.cellIDs,AllCells);
[~,matIndex_sub2] = ismember(block2_cleaned_submod.cellIDs,AllCells);

blockSubMod.TrueSynapse(1,1) = sum(sum(ConnMatrixPre_cleaned(matIndex_sub1,matIndex_sub1)));
blockSubMod.TrueSynapse(1,2) = sum(sum(ConnMatrixPre_cleaned(matIndex_sub2,matIndex_sub1)));
blockSubMod.TrueSynapse(2,2) = sum(sum(ConnMatrixPre_cleaned(matIndex_sub2,matIndex_sub2)));
blockSubMod.TrueSynapse(2,1) = sum(sum(ConnMatrixPre_cleaned(matIndex_sub1,matIndex_sub2)));
blockSubMod.Denom = [size(block1_cleaned_submod.cellIDs,2)*size(block1_cleaned_submod.cellIDs,2),...
                    size(block1_cleaned_submod.cellIDs,2)*size(block2_cleaned_submod.cellIDs,2);...
                    size(block2_cleaned_submod.cellIDs,2)*size(block1_cleaned_submod.cellIDs,2),...
                    size(block2_cleaned_submod.cellIDs,2)*size(block2_cleaned_submod.cellIDs,2)];

blockSubMod.PotentialSynapse_2(1,1) = sum(sum(potentialSynapses_2(matIndex_sub1,matIndex_sub1)));
blockSubMod.PotentialSynapse_2(1,2) = sum(sum(potentialSynapses_2(matIndex_sub2,matIndex_sub1)));
blockSubMod.PotentialSynapse_2(2,2) = sum(sum(potentialSynapses_2(matIndex_sub2,matIndex_sub2)));
blockSubMod.PotentialSynapse_2(2,1) = sum(sum(potentialSynapses_2(matIndex_sub1,matIndex_sub2)));

blockSubMod.PotentialSynapse_5(1,1) = sum(sum(potentialSynapses_5(matIndex_sub1,matIndex_sub1)));
blockSubMod.PotentialSynapse_5(1,2) = sum(sum(potentialSynapses_5(matIndex_sub2,matIndex_sub1)));
blockSubMod.PotentialSynapse_5(2,2) = sum(sum(potentialSynapses_5(matIndex_sub2,matIndex_sub2)));
blockSubMod.PotentialSynapse_5(2,1) = sum(sum(potentialSynapses_5(matIndex_sub1,matIndex_sub2)));


blockSubMod.PotentialSynapse_10(1,1) = sum(sum(potentialSynapses_10(matIndex_sub1,matIndex_sub1)));
blockSubMod.PotentialSynapse_10(1,2) = sum(sum(potentialSynapses_10(matIndex_sub2,matIndex_sub1)));
blockSubMod.PotentialSynapse_10(2,2) = sum(sum(potentialSynapses_10(matIndex_sub2,matIndex_sub2)));
blockSubMod.PotentialSynapse_10(2,1) = sum(sum(potentialSynapses_10(matIndex_sub1,matIndex_sub2)));



subplot(4,4,5)
plot([1,2,3,4],[DiagRatio(blockSubMod.TrueSynapse./blockSubMod.Denom),DiagRatio(blockSubMod.PotentialSynapse_2./blockSubMod.Denom),...
    DiagRatio(blockSubMod.PotentialSynapse_5./blockSubMod.Denom),DiagRatio(blockSubMod.PotentialSynapse_10./blockSubMod.Denom)],'-ko','LineWidth',2);
set(gca,'XTick',[1:4],'XTickLabels',[{'true'},2,5,10],'YLim',[0,6]);
box off;
daspect([1,1,1]);
%offsetAxes(gca);
xlabel('radius (\mum)');
ylabel('Sum within-module/Sum between-module')
%title('sub modules')

[~,matIndex_ABDm] = ismember([ABDr_CellIDs,ABDc_CellIDs],AllCells);
[~,matIndex_ABDi] = ismember([ABDIr_CellIDs,ABDIc_CellIDs],AllCells);


blockSubMod.ABDSynapse(1,1) = sum(sum(ConnMatrixPre_cleaned(matIndex_ABDm,matIndex_sub1)));
blockSubMod.ABDSynapse(1,2) = sum(sum(ConnMatrixPre_cleaned(matIndex_ABDm,matIndex_sub2)));
blockSubMod.ABDSynapse(2,1) = sum(sum(ConnMatrixPre_cleaned(matIndex_ABDi,matIndex_sub1)));
blockSubMod.ABDSynapse(2,2) = sum(sum(ConnMatrixPre_cleaned(matIndex_ABDi,matIndex_sub2)));
blockSubMod.ABDDenom = [size([ABDr_CellIDs,ABDc_CellIDs],2)*size(block1_cleaned_submod.cellIDs,2),...
                        size([ABDr_CellIDs,ABDc_CellIDs],2)*size(block2_cleaned_submod.cellIDs,2);...
                        size([ABDIr_CellIDs,ABDIc_CellIDs],2)*size(block1_cleaned_submod.cellIDs,2),...
                        size([ABDIr_CellIDs,ABDIc_CellIDs],2)*size(block2_cleaned_submod.cellIDs,2)];


blockSubMod.PotABDSynapse_2(1,1) = sum(sum(potentialSynapses_2(matIndex_ABDm,matIndex_sub1)));
blockSubMod.PotABDSynapse_2(1,2) = sum(sum(potentialSynapses_2(matIndex_ABDm,matIndex_sub2)));
blockSubMod.PotABDSynapse_2(2,1) = sum(sum(potentialSynapses_2(matIndex_ABDi,matIndex_sub1)));
blockSubMod.PotABDSynapse_2(2,2) = sum(sum(potentialSynapses_2(matIndex_ABDi,matIndex_sub2)));

blockSubMod.PotABDSynapse_5(1,1) = sum(sum(potentialSynapses_5(matIndex_ABDm,matIndex_sub1)));
blockSubMod.PotABDSynapse_5(1,2) = sum(sum(potentialSynapses_5(matIndex_ABDm,matIndex_sub2)));
blockSubMod.PotABDSynapse_5(2,1) = sum(sum(potentialSynapses_5(matIndex_ABDi,matIndex_sub1)));
blockSubMod.PotABDSynapse_5(2,2) = sum(sum(potentialSynapses_5(matIndex_ABDi,matIndex_sub2)));

blockSubMod.PotABDSynapse_10(1,1) = sum(sum(potentialSynapses_10(matIndex_ABDm,matIndex_sub1)));
blockSubMod.PotABDSynapse_10(1,2) = sum(sum(potentialSynapses_10(matIndex_ABDm,matIndex_sub2)));
blockSubMod.PotABDSynapse_10(2,2) = sum(sum(potentialSynapses_10(matIndex_ABDi,matIndex_sub1)));
blockSubMod.PotABDSynapse_10(2,2) = sum(sum(potentialSynapses_10(matIndex_ABDi,matIndex_sub2)));


subplot(4,4,6)
plot([1,2,3,4],[DiagRatio(blockSubMod.ABDSynapse./blockSubMod.ABDDenom,'rev'),DiagRatio(blockSubMod.PotABDSynapse_2./blockSubMod.ABDDenom,'rev'),...
    DiagRatio(blockSubMod.PotABDSynapse_5./blockSubMod.ABDDenom,'rev'),DiagRatio(blockSubMod.PotABDSynapse_10./blockSubMod.ABDDenom,'rev')],'-ko','LineWidth',2);
xlabel('radius (\mum)');
ylabel('modO-->ABD/modA -->ABD')
set(gca,'XTick',[1:4],'XTickLabels',[{'true'},2,5,10]);

box off;
daspect([1,1,1])

%%


function[n,e_index] = edgeCounts(PotentialSynapses,AllCells)
[n,e] = groupcounts(PotentialSynapses);
[~,e_index] = ismember(e,AllCells);
e_index = e_index(e_index~=0);
n = n(e_index~=0);
end

function [r] = DiagRatio(A,~)
if nargin ==1
    trueDiag = A(1,1)+A(2,2);
    offDiag = A(1,2)+A(2,1);
    r = trueDiag/offDiag;
else
    trueDiag = A(1,2)+A(2,1);
    offDiag = A(1,1)+A(2,2);
    r = trueDiag/offDiag;
end
end

function [size] = synSize(index1,index2,df_cleaned)
size = [];
for i = 1:length(index1)
    for j = 1:length(index2)
        size = [size;...
            df_cleaned.size(df_cleaned.presyn_segid == index1(i) & df_cleaned.postsyn_segid == index2(j))];
    end
end
end

% writematrix([[block1_cleaned.b1tob1size;block2_cleaned.b2tob2size],...
%     [repmat("block1",size(block1_cleaned.b1tob1size));...
%     repmat("block2",size(block2_cleaned.b2tob2size))]],'SynSize.csv')


