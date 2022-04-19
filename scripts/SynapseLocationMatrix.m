% synapse locations based on the modular connectome


load MatOrder_CO_top500_2blocks_gamma038_08062020.mat
load cellIDType_CO_top500_2blocks_gamma038_08062020.mat
load AllCells.mat
load ConnMatrixPre_cleaned.mat
load df_cleaned.mat


cellIDType = cellstr(cellIDType_CO_top500_2blocks_gamma038_08062020);

block1_cleaned.cellIDs = MatOrder_CO_top500_2blocks_gamma038_08062020(strcmp(cellIDType,'_Axl_'));
block2_cleaned.cellIDs = [MatOrder_CO_top500_2blocks_gamma038_08062020(strcmp(cellIDType,'_Int_')),...
                          MatOrder_CO_top500_2blocks_gamma038_08062020(strcmp(cellIDType,'_DOs_'))];
vspnCells = MatOrder_CO_top500_2blocks_gamma038_08062020(strcmp(cellIDType,'vSPNs'));
ABDm_cells = MatOrder_CO_top500_2blocks_gamma038_08062020(strcmp(cellIDType,'ABD_m'));
ABDi_cells = MatOrder_CO_top500_2blocks_gamma038_08062020(strcmp(cellIDType,'ABD_i'));

block1_cleaned.color = [1,0.5,0];
block2_cleaned.color = [0,0.5,1];


OrderedCells = [block1_cleaned.cellIDs,block2_cleaned.cellIDs,vspnCells,...
                ABDm_cells,ABDi_cells];
            
[~,OrderedIndex] =ismember(OrderedCells,AllCells);
connMat = ConnMatrixPre_cleaned(OrderedIndex,OrderedIndex);


%%
meanLengthMat = zeros(size(OrderedCells,2));
meanLengthMat_normalized = zeros(size(OrderedCells,2));

for i = 1:size(OrderedCells,2)
    if isExistReRoot(OrderedCells(i)) ~= 0
        [temp1,temp2] = SynapticPartners(OrderedCells(i),1,df_cleaned);
        
        temp1 = temp1(temp1<1e5);
        partnerIndex = ismember(temp1,OrderedCells);
        partnerCoord = PrePartnerCoordinates(temp2(partnerIndex),df_cleaned);
        if ~isempty(partnerCoord)
            partnerID{i} = temp1(partnerIndex);
            partnerPathlength{i} = PathLengthToCoordinate(partnerCoord,OrderedCells(i));
            
            lengthMat =  nan(22,size(OrderedCells,2));
            colInd = [];
            for j = 1:length(partnerID{i})
                colIndex = find(OrderedCells == partnerID{i}(j,1));
                colInd = [colInd,colIndex];
                [a,b] = histcounts(colInd,1:length(OrderedCells));
                rowIndex = a(b==colIndex);
                lengthMat(rowIndex+1,colIndex) = partnerPathlength{i}(j,1);
            end
            meanLengthMat(i,:) = nanmean(lengthMat,1);
            meanLengthMat_normalized(i,:) = nanmean(lengthMat,1)./TreeLength(OrderedCells(i));
        else
            continue
        end
        
        clear temp1
        clear temp2
        clear partnerIndex
        clear partnerCoord
    end
end

%%
figure;
subplot(1,2,1)
cspy(meanLengthMat_normalized,'Colormap',colorcet('R3','N',5),'Levels',5,'MarkerSize',8);
colorbar();
line([length(block1_cleaned.cellIDs),length(block1_cleaned.cellIDs)],[0,size(meanLengthMat_normalized,1)],'color','k');
hold on;
line([length(block1_cleaned.cellIDs),length(block1_cleaned.cellIDs)]+length(block2_cleaned.cellIDs),[0,size(meanLengthMat_normalized,1)],'color','k');
line([length(block1_cleaned.cellIDs),length(block1_cleaned.cellIDs)]+length(block2_cleaned.cellIDs)-34,[0,size(meanLengthMat_normalized,1)],'color','k');
line([length(block1_cleaned.cellIDs),length(block1_cleaned.cellIDs)]+length(block2_cleaned.cellIDs)+ length(vspnCells)...
    ,[0,size(meanLengthMat_normalized,1)],'color','k');

line([0,size(meanLengthMat_normalized,1)],[length(block1_cleaned.cellIDs),length(block1_cleaned.cellIDs)],'color','k');
line([0,size(meanLengthMat_normalized,1)],[length(block1_cleaned.cellIDs),length(block1_cleaned.cellIDs)]+length(block2_cleaned.cellIDs),'color','k');
line([0,size(meanLengthMat_normalized,1)],[length(block1_cleaned.cellIDs),length(block1_cleaned.cellIDs)]+length(block2_cleaned.cellIDs)-34,'color','k');
line([0,size(meanLengthMat_normalized,1)],...
    [length(block1_cleaned.cellIDs),length(block1_cleaned.cellIDs)]+length(block2_cleaned.cellIDs)+ length(vspnCells),'color','k');
box on;
title('Normalized pathlength to avg synapse');

subplot(1,2,2)
cspy(connMat,'Colormap',colorcet('R3','N',5),'Levels',5,'MarkerSize',8);
line([length(block1_cleaned.cellIDs),length(block1_cleaned.cellIDs)],[0,size(meanLengthMat_normalized,1)],'color','k');
hold on;
line([length(block1_cleaned.cellIDs),length(block1_cleaned.cellIDs)]+length(block2_cleaned.cellIDs),[0,size(meanLengthMat_normalized,1)],'color','k');
line([length(block1_cleaned.cellIDs),length(block1_cleaned.cellIDs)]+length(block2_cleaned.cellIDs)-34,[0,size(meanLengthMat_normalized,1)],'color','k');
line([length(block1_cleaned.cellIDs),length(block1_cleaned.cellIDs)]+length(block2_cleaned.cellIDs)+ length(vspnCells)...
    ,[0,size(meanLengthMat_normalized,1)],'color','k');

line([0,size(meanLengthMat_normalized,1)],[length(block1_cleaned.cellIDs),length(block1_cleaned.cellIDs)],'color','k');
line([0,size(meanLengthMat_normalized,1)],[length(block1_cleaned.cellIDs),length(block1_cleaned.cellIDs)]+length(block2_cleaned.cellIDs)-34,'color','k');
line([0,size(meanLengthMat_normalized,1)],[length(block1_cleaned.cellIDs),length(block1_cleaned.cellIDs)]+length(block2_cleaned.cellIDs),'color','k');
line([0,size(meanLengthMat_normalized,1)],...
    [length(block1_cleaned.cellIDs),length(block1_cleaned.cellIDs)]+length(block2_cleaned.cellIDs)+ length(vspnCells),'color','k');
box on;
title('No. synapses');

%% arrange block2 by RC position

OM.cellIDs = block2_cleaned.cellIDs;
OM.allCellsIDs = [block2_cleaned.cellIDs,vspnCells,ABDm_cells,ABDi_cells];
OM.origins = getOrigin(OM.cellIDs);
[~,OM.sortedOriginIndex] = sort(OM.origins(:,2));

OM.meanLengthMat = meanLengthMat_normalized(252:end,252:end);
OM.connMat = connMat(252:end,252:end);

% block2 cell IDs 252:289

sortedIndex = [OM.sortedOriginIndex',290:358];
OM.sortedMeanMat = OM.meanLengthMat(sortedIndex,sortedIndex);
OM.sortedConnMat = OM.connMat(sortedIndex,sortedIndex);
%%

figure;
subplot(1,2,1)
cspy(OM.sortedMeanMat,'Colormap',colorcet('R3','N',5),'Levels',5,'MarkerSize',8);
box on;
line([length(block2_cleaned.cellIDs),length(block2_cleaned.cellIDs)],[0,size(OM.sortedMeanMat,1)],'color','k');
line([length(block2_cleaned.cellIDs),length(block2_cleaned.cellIDs)]-34,[0,size(OM.sortedMeanMat,1)],'color','k');
line([length(block2_cleaned.cellIDs),length(block2_cleaned.cellIDs)]+length(vspnCells),[0,size(OM.sortedMeanMat,1)],'color','k');

line([0,size(OM.sortedMeanMat,1)],[length(block2_cleaned.cellIDs),length(block2_cleaned.cellIDs)],'color','k');
line([0,size(OM.sortedMeanMat,1)],[length(block2_cleaned.cellIDs),length(block2_cleaned.cellIDs)]-34,'color','k');
line([0,size(OM.sortedMeanMat,1)],[length(block2_cleaned.cellIDs),length(block2_cleaned.cellIDs)]+length(vspnCells),'color','k');
colorbar()
title({'modO only','Norm. pathlength to avg synapse - RC sorted'});

subplot(1,2,2)
cspy(OM.sortedConnMat,'Colormap',colorcet('R3','N',5),'Levels',5,'MarkerSize',8);
box on;
line([length(block2_cleaned.cellIDs),length(block2_cleaned.cellIDs)],[0,size(OM.sortedMeanMat,1)],'color','k');
line([length(block2_cleaned.cellIDs),length(block2_cleaned.cellIDs)]-34,[0,size(OM.sortedMeanMat,1)],'color','k');
line([length(block2_cleaned.cellIDs),length(block2_cleaned.cellIDs)]+length(vspnCells),[0,size(OM.sortedMeanMat,1)],'color','k');

line([0,size(OM.sortedMeanMat,1)],[length(block2_cleaned.cellIDs),length(block2_cleaned.cellIDs)],'color','k');
line([0,size(OM.sortedMeanMat,1)],[length(block2_cleaned.cellIDs),length(block2_cleaned.cellIDs)]-34,'color','k');
line([0,size(OM.sortedMeanMat,1)],[length(block2_cleaned.cellIDs),length(block2_cleaned.cellIDs)]+length(vspnCells),'color','k');
title({'modO - only','No.of synapses - RC sorted'});

%%

OM.sortedMeanMat_T = OM.sortedMeanMat;
OM.sortedMeanMat_T(OM.sortedMeanMat>0.4) = 0.4;
figure;
subplot(1,2,1)
cspy(OM.sortedMeanMat_T,'Colormap',colorcet('R3','N',4),'Levels',4,'MarkerSize',12);
box on;
line([length(block2_cleaned.cellIDs),length(block2_cleaned.cellIDs)],[0,size(OM.sortedMeanMat,1)],'color','k');
line([length(block2_cleaned.cellIDs),length(block2_cleaned.cellIDs)]-34,[0,size(OM.sortedMeanMat,1)],'color','k');
line([length(block2_cleaned.cellIDs),length(block2_cleaned.cellIDs)]+length(vspnCells),[0,size(OM.sortedMeanMat,1)],'color','k');

line([0,size(OM.sortedMeanMat,1)],[length(block2_cleaned.cellIDs),length(block2_cleaned.cellIDs)],'color','k');
line([0,size(OM.sortedMeanMat,1)],[length(block2_cleaned.cellIDs),length(block2_cleaned.cellIDs)]-34,'color','k');
line([0,size(OM.sortedMeanMat,1)],[length(block2_cleaned.cellIDs),length(block2_cleaned.cellIDs)]+length(vspnCells),'color','k');
colorbar()
title({'modO only','Norm. pathlength to avg synapse - RC sorted'});
caxis([0,0.4]);

OM.sortedConnMat_T = OM.sortedConnMat;
OM.sortedConnMat_T(OM.sortedConnMat>5) = 5;
subplot(1,2,2)
cspy(OM.sortedConnMat_T,'Colormap',colorcet('R3','N',5),'Levels',5,'MarkerSize',12);
box on;
line([length(block2_cleaned.cellIDs),length(block2_cleaned.cellIDs)],[0,size(OM.sortedMeanMat,1)],'color','k');
line([length(block2_cleaned.cellIDs),length(block2_cleaned.cellIDs)]-34,[0,size(OM.sortedMeanMat,1)],'color','k');
line([length(block2_cleaned.cellIDs),length(block2_cleaned.cellIDs)]+length(vspnCells),[0,size(OM.sortedMeanMat,1)],'color','k');

line([0,size(OM.sortedMeanMat,1)],[length(block2_cleaned.cellIDs),length(block2_cleaned.cellIDs)],'color','k');
line([0,size(OM.sortedMeanMat,1)],[length(block2_cleaned.cellIDs),length(block2_cleaned.cellIDs)]-34,'color','k');
line([0,size(OM.sortedMeanMat,1)],[length(block2_cleaned.cellIDs),length(block2_cleaned.cellIDs)]+length(vspnCells),'color','k');
title({'modO - only','No.of synapses - RC sorted'});

