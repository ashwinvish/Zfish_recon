% plot K values from python


clc;
clear;

% vpni = load('fitCellsBasedOnLocationOthers.mat');
% vpni.zCoord(:,1) = 0.798005* vpni.zCoord(:,1); % RC
% vpni.zCoord(:,2) = 0.798005* vpni.zCoord(:,2); % ML
% vpni.zCoord(:,3) = 2*0.798005* vpni.zCoord(:,3); % DV
% vpni.rhombomere = isRhombomere('',[vpni.zCoord(:,2),vpni.zCoord(:,1),vpni.zCoord(:,3)]);
% save('fitCellsBasedOnLocationOthers_withRhombomere.mat','vpni');


% load only raw K values
S = load('/Users/ashwin/Google Drive/Zfish/ZFishPaper-2/forNature/Nature-FigureList/newMatfiles/slopes_CO_top500_08062020.mat');
load cellIDType_CO_top500_2blocks_gamma038_08062020.mat
load MatOrder_CO_top500_2blocks_gamma038_08062020.mat


% load 1000 variants of K values

%S = load('/Users/ashwin/Google Drive/Zfish/ZFishPaper-2/forNature/Nature-FigureList/newMatFiles/slopes_many_CO_top500_08062020.mat');


Int_rm = S.Int;
Vest_rm = S.Vest;
ABDm_rm = S.ABDm;
ABDi_rm = S.ABDi;

cellRhombomeres = isRhombomere(MatOrder_CO_top500_2blocks_gamma038_08062020);

r46Int = [S.Kvalues(:,logical(cellRhombomeres.r4)),S.Kvalues(:,logical(cellRhombomeres.r5)),S.Kvalues(:,logical(cellRhombomeres.r6))];
r46Int = r46Int(r46Int>0);
r78Int = S.Kvalues(:,logical(cellRhombomeres.r7));
r78Int = r78Int(r78Int>0);

% calculate means
model_mean = [mean(Int_rm),mean(Vest_rm), mean(ABDm_rm),mean(ABDi_rm)];
model_SD = [std(Int_rm)/sqrt(size(Int_rm,1)),std(Vest_rm)/sqrt(size(Vest_rm,1)) std(ABDm_rm)/sqrt(size(ABDm_rm,1)),std(ABDi_rm)/sqrt(size(ABDi_rm,1))];

macaque_mean = [2.74,NaN,6.2,4.6];
macaque_SD = [0.17,NaN,0.34,0.25];

cat_mean = [NaN,NaN, 5.61,  12.41];
cat_SD = [NaN, NaN,0.15 ,0.92];

GF_mean = [2.8,0.17,6.76,8.37];
GF_SD = [0.57 ,0.02, 0.54,0.39];


% experimental data
load 'fitCellsBasedOnLocationOthers_withRhombomere.mat' % vpni
abd = load('fitCellsBasedOnLocationAbd');
abdi = load('fitCellsBasedOnLocationAbdi');
do = load('fitCellsBasedOnLocationDO');


% experimental cutoffs
varCut = 1;
r2CutVPNI = 0.4;
r2CutABD = 0.4;
nMax = size(vpni.slopes(:,1),1);

exptSlopes = [vpni.slopes(:,1),...
    [do.slopes;NaN(nMax-size(do.slopes,1),1)],...
    [abd.slopes(:,1);NaN(nMax-size(abd.slopes,1),1)],...
    [abdi.slopes(:,2);NaN(nMax-size(abdi.slopes,1),1)]];

exptSlopes(~(vpni.var2explain(:,1)>varCut & vpni.modelStats(:,1)>r2CutVPNI),1) = NaN;
exptSlopes(do.var2explain<varCut,2)=NaN;
exptSlopes(~(abd.var2explain(:,1)>varCut & abd.modelStats(:,1)>r2CutABD),3) = NaN;
exptSlopes(~(abdi.var2explain(:,1)>varCut & abdi.modelStats(:,2)>r2CutABD),4) = NaN;

exptVPNI = exptSlopes(:,1);
exptDO = removeOutliers(exptSlopes(:,2));
exptABD = removeOutliers(exptSlopes(:,3));
exptABDi  = removeOutliers(exptSlopes(:,4));

exptVPNIr46 = [exptVPNI(logical(vpni.rhombomere.r4));...
    exptVPNI(logical(vpni.rhombomere.r5));...
    exptVPNI(logical(vpni.rhombomere.r6))];
exptVPNIr46 = removeOutliers(exptVPNIr46);

exptVPNIr78 = exptVPNI(logical(vpni.rhombomere.r7));
exptVPNIr78 = removeOutliers(exptVPNIr78);

exptVPNI = removeOutliers( exptSlopes(:,1));

% scale to Goldfish
sf = 2.8 / nanmean(exptVPNIr78);
sf = 1;

%% 
figure; 


modelData = NaN(size(r46Int,2),5);
modelData(1:length(r46Int),1) = r46Int';
modelData(1:length(r78Int),2) = r78Int';
modelData(1:length(Vest_rm),3) = Vest_rm';
modelData(1:length(ABDm_rm),4) = ABDm_rm';
modelData(1:length(ABDi_rm),5) = ABDi_rm';

exptData = NaN(sum(~isnan(exptVPNIr46)),5);
exptData(1:sum(~isnan(exptVPNIr46)),1) = exptVPNIr46(~isnan(exptVPNIr46));
exptData(1:sum(~isnan(exptVPNIr78)),2) = exptVPNIr78(~isnan(exptVPNIr78));
exptData(1:sum(~isnan(exptDO)),3) = exptDO(~isnan(exptDO));
exptData(1:sum(~isnan(exptABD)),4) = exptABD(~isnan(exptABD));
exptData(1:sum(~isnan(exptABDi)),5) = exptABDi(~isnan(exptABDi));


violinPlot(modelData,'histOri','left','color','k','widthDiv',[2,1],'showMM',1,'histOpt',0);
violinPlot(sf*exptData,'histOri','right','color','r','widthDiv',[2,2],'showMM',1,'histOpt',0);
set(gca,'YGrid','on','YMinorGrid','on','XTick',[1,2,3,4,5],'XTickLabels',{'r46','r78','Vest','ABDm','ABDi'})



%%
figure;
subplot(4,4,1)

violin(1,Int_rm,'LineWidth',2,'style',1,'FaceColor',[0,0.5,1],'side','both');
hold on
violin(2,Vest_rm,'LineWidth',2,'style',1,'FaceColor','r','side','both');
violin(3,ABDm_rm,'LineWidth',2,'style',1,'FaceColor',[0.0196,0.6000,0.2824],'side','both');
violin(4,ABDi_rm,'LineWidth',2,'style',1,'FaceColor',[0.9216,0,0.5451],'side','both');

plot(1,nanmean(Int_rm),'k*')
plot(2,nanmean(Vest_rm),'k*')
plot(3,nanmean(ABDm_rm),'k*')
plot(4,nanmean(ABDi_rm),'k*')

set(gca,'YGrid','on','YMinorGrid','on','XTick',[1,2,3,4],'XTickLabels',{'Integ','Vest','ABDm','ABDi'})
title('Model (rawK)')


subplot(4,4,2)

violin(1,sf*exptVPNI(~isnan(exptVPNI)),'LineWidth',2,'style',1,'FaceColor',[0,0.5,1],'side','both');
hold on
violin(2,sf*exptDO(~isnan(exptDO)),'LineWidth',2,'style',1,'FaceColor','r','side','both','scaling',0.05);
violin(3,sf*exptABD(~isnan(exptABD)),'LineWidth',2,'style',1,'FaceColor',[0.0196,0.6000,0.2824],'side','both');
violin(4,sf*exptABDi(~isnan(exptABDi)),'LineWidth',2,'style',1,'FaceColor',[0.9216,0,0.5451],'side','both');

plot(1,nanmean(sf *exptVPNI),'ro')
plot(2,nanmean(sf *exptDO),'ro')
plot(3,nanmean(sf *exptABD),'ro')
plot(4,nanmean(sf *exptABDi),'ro')

set(gca,'YGrid','on','YMinorGrid','on','XTick',[1,2,3,4],'XTickLabels',{'Integ','Vest','ABDm','ABDi'})
title('Experimental')


subplot(4,4,3)

plot([1,2,3,4],[nanmean(Int_rm),nanmean(Vest_rm),nanmean(ABDm_rm),nanmean(ABDi_rm)],'-k*')
hold on;
plot([1,2,3,4],[nanmean(sf *exptVPNI),nanmean(sf *exptDO),nanmean(sf *exptABD),nanmean(sf *exptABDi)],'-ro')
set(gca,'YGrid','on','YMinorGrid','on','XTick',[1,2,3,4],'XTickLabels',{'Integ','Vest','ABDm','ABDi'})
legend({'model','expt'});
box off;

subplot(4,4,4)

violin(1,Int_rm,'LineWidth',2,'style',1,'FaceColor',[0,0.5,1],'side','both');
hold on
violin(2,Vest_rm,'LineWidth',2,'style',1,'FaceColor','r','side','both');
violin(3,ABDm_rm,'LineWidth',2,'style',1,'FaceColor',[0.0196,0.6000,0.2824],'side','both');
violin(4,ABDi_rm,'LineWidth',2,'style',1,'FaceColor',[0.9216,0,0.5451],'side','both');

plot(1,nanmean(Int_rm),'k*')
plot(2,nanmean(Vest_rm),'k*')
plot(3,nanmean(ABDm_rm),'k*')
plot(4,nanmean(ABDi_rm),'k*')

plot(1,nanmean(sf *exptVPNI),'ro')
plot(2,nanmean(sf *exptDO),'ro')
plot(3,nanmean(sf *exptABD),'ro')
plot(4,nanmean(sf *exptABDi),'ro')

set(gca,'YGrid','on','YMinorGrid','on','XTick',[1,2,3,4],'XTickLabels',{'Integ','Vest','ABDm','ABDi'})
box off;
title('Model Vs. Expt');

subplot(4,4,5)
violin(1,r46Int,'LineWidth',2,'style',1,'FaceColor',[0,0.5,1],'side','both');
hold on;
violin(2,r78Int,'LineWidth',2,'style',1,'FaceColor',[0,0.5,1],'side','both');
violin(3,Vest_rm,'LineWidth',2,'style',1,'FaceColor','r','side','both');
violin(4,ABDm_rm,'LineWidth',2,'style',1,'FaceColor',[0.0196,0.6000,0.2824],'side','both');
violin(5,ABDi_rm,'LineWidth',2,'style',1,'FaceColor',[0.9216,0,0.5451],'side','both');

plot(1,nanmean(r46Int),'k*')
plot(2,nanmean(r78Int),'k*')
plot(3,nanmean(Vest_rm),'k*')
plot(4,nanmean(ABDm_rm),'k*')
plot(5,nanmean(ABDi_rm),'k*')

set(gca,'YGrid','on','YMinorGrid','on','XTick',[1,2,3,4,5],'XTickLabels',{'r46','r78','Vest','ABDm','ABDi'})
title('Model (rawK)')

subplot(4,4,6)

violin(1,sf*exptVPNIr46(~isnan(exptVPNIr46)),'LineWidth',2,'style',1,'FaceColor',[0,0.5,1],'side','both');
hold on
violin(2,sf*exptVPNIr78(~isnan(exptVPNIr78)),'LineWidth',2,'style',1,'FaceColor',[0,0.5,1],'side','both');
violin(3,sf*exptDO(~isnan(exptDO)),'LineWidth',2,'style',1,'FaceColor','r','side','both','scaling',0.05);
violin(4,sf*exptABD(~isnan(exptABD)),'LineWidth',2,'style',1,'FaceColor',[0.0196,0.6000,0.2824],'side','both');
violin(5,sf*exptABDi(~isnan(exptABDi)),'LineWidth',2,'style',1,'FaceColor',[0.9216,0,0.5451],'side','both');

plot(1,nanmean(sf *exptVPNIr46),'ro')
plot(2,nanmean(sf *exptVPNIr78),'ro')
plot(3,nanmean(sf *exptDO),'ro')
plot(4,nanmean(sf *exptABD),'ro')
plot(5,nanmean(sf *exptABDi),'ro')

set(gca,'YGrid','on','YMinorGrid','on','XTick',[1,2,3,4,5],'XTickLabels',{'r46','r78','Vest','ABDm','ABDi'})
title('Experimental')

subplot(4,4,7)
plot([1,2,3,4,5],[nanmean(r46Int),nanmean(r78Int),nanmean(Vest_rm),nanmean(ABDm_rm),nanmean(ABDi_rm)],'-k*');
hold on;
plot([1,2,3,4,5],[nanmean(sf *exptVPNIr46),nanmean(sf *exptVPNIr78),nanmean(sf *exptDO),nanmean(sf *exptABD),nanmean(sf *exptABDi)],'-ro');
set(gca,'YGrid','on','YMinorGrid','on','XTick',[1,2,3,4,5],'XTickLabels',{'r46','r78','Vest','ABDm','ABDi'})
legend({'model','expt'});
box off;

%%
subplot(4,4,3)

errorbar(model_mean,model_SD,'-ko');
hold on;
errorbar(macaque_mean,macaque_SD,'-.o','LineWidth',2);
errorbar(cat_mean,cat_SD,'-.o','LineWidth',2);
errorbar(GF_mean,GF_SD,'-.o','LineWidth',2);
box off;
axis square;
set(gca,'XTick',[1,2,3,4],'XTickLabels',{'Int','vest','ABDm','ABDi'},'YLim',[0,15]);
ylabel('K')
offsetAxes(gca);
legend({'model','mac','cat','gf'},'Location','bestoutside');



%% find top50 intervals

load cellIDType_CO_top500_2blocks_gamma038_08062020.mat
load MatOrder_CO_top500_2blocks_gamma038_08062020.mat
load cellIDType_CO_top500_2blocks_gamma038_08062020.mat

cellType = cellstr(cellIDType_CO_top500_2blocks_gamma038_08062020);
IntCells = strcmp(cellType,'_Int_');
VestCells = strcmp(cellType,'_DOs_');
AxlCells = strcmp(cellType,'_Axl_');

OMcells = logical(IntCells+VestCells);
OMCellIDs = MatOrder_CO_top500_2blocks_gamma038_08062020(OMcells);

fileName = 'ranking_noMirror_cleaned_squared.csv';
fid = fopen(fileName);
% cellID,fanOut,fanIn,toIntegrator,fromInteg,toMotor,class,centrality
SpectralWithRS.complete = textscan(fid,'%d %d %d %d %d %d %d %s %f','Delimiter',',');
fclose(fid);
ranksToConsider = 1:sum(SpectralWithRS.complete{9}>0);
cellIDsinECorder = SpectralWithRS.complete{2}(ranksToConsider);

[lia,lob] = ismember(OMCellIDs,cellIDsinECorder);




%top50intervals = [164,181,204,230,257,287,319,349,378,403,434,458]; %458


% remove outliers <0.3 percentile and >99.7percentile
function[population_rm] = removeOutliers(population)
population = population(:);
low = prctile(population,5);
high = prctile(population,95);
population_rm = population((population>low) & (population<high));
end

