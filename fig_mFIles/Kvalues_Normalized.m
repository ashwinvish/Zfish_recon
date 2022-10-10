clc;
clear;

load ('../matFiles/cellIDType_CO_top500_2blocks_gamma038_08062020.mat');
load('../matFiles/MatOrder_CO_top500_2blocks_gamma038_08062020.mat');
cellIDType = cellstr(cellIDType_CO_top500_2blocks_gamma038_08062020);


S = load('../matFiles/slopes_many_CO_top500_unscaled_08062020.mat');

Int_rm = removeOutliers(S.Int,[1,99]);
Vest_rm = removeOutliers(S.Vest,[1,99]);
ABDm_rm = removeOutliers(S.ABDm,[1,99]);
ABDi_rm = removeOutliers(S.ABDi,[1,99]);

Int_rm(Int_rm==0) = NaN;
Vest_rm(Vest_rm==0) = NaN;
ABDm_rm(ABDm_rm==0) = NaN;
ABDi_rm (ABDi_rm==0) = NaN;

cellRhombomeres = isRhombomere(MatOrder_CO_top500_2blocks_gamma038_08062020(strcmp(cellIDType,'_Int_')));

r46Int = [S.Int(:,logical(cellRhombomeres.r4)),S.Int(:,logical(cellRhombomeres.r5)),S.Int(:,logical(cellRhombomeres.r6))];
r78Int = S.Int(:,logical(cellRhombomeres.r7));


% calculate means
macaque_mean(1) = (51*3.2 + 19*2.4 + 29*2.5)/(51+19+29);
macaque_SE(1) = sqrt(1.3^2/51+1.6^2/19+1.0^2/29)/3;

macaque_mean(2) = NaN;
macaque_SE(2) = NaN;

macaque_mean(3) = (81*6.2 + 46*5.9 + 22*3.4)/(81+46+22);
macaque_SE(3) = sqrt(3.09^2/81 + 1.65^2/46 + 0.87^2/22)/3;

macaque_mean(4) = (36*4.6 + 24*5.3 + 18*4.92)/ (36+24+18);
macaque_SE(4) =  sqrt(1.5^2/36 + 1.3^2/24 + 1.74^2/18)/3;


cat_mean(1) = (6*7.3 + 10*8.2 + 64*3.7)/(6+10+64);
cat_SE(1) = sqrt(2.7^2/6 + 2.9^2/10 + 3.13^2/64) / 3;

cat_mean(2) = NaN;
cat_SE(2) = NaN;

cat_mean(3) = ( 40*8.7 + 104*6.7 + 51*3.13 + 64* 8.8 + 118*6.25 ) / (40+104+51+64+118);
cat_SE(3) = sqrt(2.5^2/40 + 0.19^2 + 0.15^2 + 2.2^2/64 + 3/05^2/118) / 5;

cat_mean(4) = (43*12.01 + 24*7.01 + 50*7) / (43+24+50);
cat_SE(4) = sqrt(3.1^2/43 + 1.64^2/24 + 2.69^2/50)/3;

GF_mean(1) = (44*2.8 + 11*1.67)/(44+11);
GF_SE(1) = 2.8*.205 + 0.51/sqrt(11);

GF_mean(2) = 0.1;
GF_SE(2) = 0.02;

GF_mean(3) = (7.13*12 + 6.32*10) / 22;
GF_SE(3) = sqrt(3.4^2/12 + 1.4^2/10)/2;

GF_mean(4) = 8.37;
GF_SE(4) = 1.71 / sqrt(19);

macaque_SF = GF_mean(1)/macaque_mean(1);
cat_SF = GF_mean(1)/cat_mean(1);

% 1000 models
modelData_1000 = NaN(size(Int_rm,1)*size(Int_rm,2),5);
modelData_1000(1:length(r46Int(:)),1) = r46Int(:);
modelData_1000(1:length(r78Int(:)),2) = r78Int(:);
modelData_1000(1:length(Vest_rm(:)),3) = Vest_rm(:);
modelData_1000(1:length(ABDm_rm(:)),4) = ABDm_rm(:);
modelData_1000(1:length(ABDi_rm(:)),5) = ABDi_rm(:);


% experimental
load ('../matFiles/slopesThresh_09102020.csv'); % 4 - ABD; 5 - ABDi, 6 - Do, 7 - threshold; 8 - slopes

exptABD = slopesThresh_09102020(logical(slopesThresh_09102020(:,4)),8);
exptABDi = slopesThresh_09102020(logical(slopesThresh_09102020(:,5)),8);
exptDO = slopesThresh_09102020(logical(slopesThresh_09102020(:,6)),8);
exptDO = load('../matFiles/fitCellsBasedOnLocationDO.mat');

exptDO.slopes(exptDO.var2explain<1)=NaN; % varCut = 1;

exptVPNI_index = 1-(slopesThresh_09102020(:,4)+slopesThresh_09102020(:,5)+slopesThresh_09102020(:,6));
exptVPNI = slopesThresh_09102020(logical(exptVPNI_index),8);

exptVPNIcoords = [slopesThresh_09102020(logical(exptVPNI_index),2),...
    slopesThresh_09102020(logical(exptVPNI_index),1),...
    slopesThresh_09102020(logical(exptVPNI_index),3)];

exptVPNIRhombomere = isRhombomere('',exptVPNIcoords);
exptr46 = exptVPNI(logical(exptVPNIRhombomere.r4+exptVPNIRhombomere.r5+exptVPNIRhombomere.r6));
exptr78 = exptVPNI(logical(exptVPNIRhombomere.r7));

exptData = NaN(size(exptVPNI,1),5);

exptData(1:length(exptr46),1) = exptr46;
exptData(1:length(exptr78),2) = exptr78;
exptData(1:length(exptDO.slopes),3) = exptDO.slopes;
exptData(1:length(exptABD),4) = exptABD;
exptData(1:length(exptABDi),5) = exptABDi;

expt_SF = GF_mean(1)/ nanmean([exptData(:,1);exptData(:,2)]);
model_SF = GF_mean(1)/ nanmean([modelData_1000(:,1);modelData_1000(:,2)]);


exptData_1_99 = removeOutliers(exptData,[1,99]);
exptData_1_99(exptData_1_99==0) = NaN;

expt_SF_1_99 = GF_mean(1)/ nanmean([exptData_1_99(:,1);exptData_1_99(:,2)]);

% raw connectome
S_raw_actual = load('../matFiles/slopes_CO_top500_unscaled_08062020.mat');

S_raw_actual.r46 = [S_raw_actual.Int(logical(cellRhombomeres.r4)),...
    S_raw_actual.Int(logical(cellRhombomeres.r5)),...
    S_raw_actual.Int(logical(cellRhombomeres.r6))];
S_raw_actual.r78 = S_raw_actual.Int(logical(cellRhombomeres.r7));

S_raw_actual.sf = GF_mean(1)/nanmean([S_raw_actual.r46,S_raw_actual.r78]);

rawData = NaN(size(S_raw_actual.r46,2),5);

rawData(1:length(S_raw_actual.r46),1) = S_raw_actual.r46(:);
rawData(1:length(S_raw_actual.r78),2) = S_raw_actual.r78(:);
rawData(1:length(S_raw_actual.Vest),3) = S_raw_actual.Vest(:);
rawData(1:length(S_raw_actual.ABDm),4) = S_raw_actual.ABDm(:);
rawData(1:length(S_raw_actual.ABDi),5) = S_raw_actual.ABDi(:);
%%

cols = cbrewer('qual','Dark2',5);

subplot(1,2,1)

h = boxplot( S_raw_actual.sf*rawData(:),...
    [ones(numel(rawData(:,1:2)),1);2*ones(numel(rawData(:,3)),1);3*ones(numel(rawData(:,4)),1);4*ones(numel(rawData(:,5)),1)],...
    'Positions',[1,2,3,4],'colors','k','Jitter',0.1,'symbol','o');
set(h,{'linew'},{2},'lineStyle','-')
hold on;

errorbar([0.95,1.95,2.95,3.95],[mean([S_raw_actual.sf*S_raw_actual.r46,S_raw_actual.sf*S_raw_actual.r78]),mean(S_raw_actual.sf*S_raw_actual.Vest),mean(S_raw_actual.sf*S_raw_actual.ABDm),mean(S_raw_actual.sf*S_raw_actual.ABDi)],...
    [std(S_raw_actual.sf*S_raw_actual.Int)./sqrt(length(S_raw_actual.Int)),std(S_raw_actual.sf*S_raw_actual.Vest)./sqrt(length(S_raw_actual.Vest)),std(S_raw_actual.sf*S_raw_actual.ABDm)./sqrt(length(S_raw_actual.ABDm)),std(S_raw_actual.sf*S_raw_actual.ABDi)./sqrt(length(S_raw_actual.ABDi))],...
    'k^','LineWidth',2,'CapSize',0,'MarkerFaceColor','w','MarkerSize',9);

errorbar([1.05,2.05,3.05,4.05],macaque_SF*macaque_mean,macaque_SF*macaque_SE,'^','color',cols(1,:),'LineWidth',2,'CapSize',0,'MarkerFaceColor','w','MarkerSize',9);
errorbar([0.85,1.85,2.85,3.85],cat_SF*cat_mean,cat_SF*cat_SE,'^','color',cols(2,:),'LineWidth',2,'CapSize',0,'MarkerFaceColor','w','MarkerSize',9);
errorbar([1.1,2.1,3.1,4.1],GF_mean,GF_SE,'^','color',cols(3,:),'LineWidth',2,'CapSize',0,'MarkerFaceColor','w','MarkerSize',9);

errorbar([0.95,2,3.5],[nanmean([expt_SF_1_99*exptData_1_99(:,1);expt_SF_1_99*exptData_1_99(:,2)]),...
    nanmean(expt_SF_1_99*exptData_1_99(:,3)),...
    nanmean([expt_SF_1_99*exptData_1_99(:,4);expt_SF_1_99*exptData_1_99(:,5)])],...
    [nanstd([expt_SF_1_99*exptData_1_99(:,1);expt_SF_1_99*exptData_1_99(:,2)])./sqrt(sum(~isnan([expt_SF_1_99*exptData_1_99(:,1);expt_SF_1_99*exptData_1_99(:,2)]))),...
    nanstd(expt_SF_1_99*exptData_1_99(:,3))./sqrt(sum(~isnan(exptData_1_99(:,3)))),...
    nanstd([expt_SF_1_99*exptData_1_99(:,4);expt_SF_1_99*exptData_1_99(:,5)])./...
    (sqrt(sum(~isnan(exptData_1_99(:,4))) + sum(~isnan(exptData_1_99(:,5)))))],...
    '^','color',cols(4,:),'lineWidth',2,'CapSize',0,'MarkerFaceColor','w','MarkerSize',9);


daspect([1,2,1]);
box off;
set(gca,'XTickLabels',{'VPNI','Vest/DO','ABDm','ABDi'});
ylabel('Norm. eye position sensitivity (k)')
legend({'Model','Macaque','Cat','GF','Expt.'})
% errorbar([mean([S_raw.sf*S_raw.r46,S_raw.sf*S_raw.r78]),mean(S_raw.sf*S_raw.Vest),mean(S_raw.sf*S_raw.ABDm),mean(S_raw.sf*S_raw.ABDi)],...
%     [std(S_raw.sf*S_raw.Int)./sqrt(length(S_raw.Int)),std(S_raw.sf*S_raw.Vest)./sqrt(length(S_raw.Vest)),std(S_raw.sf*S_raw.ABDm)./sqrt(length(S_raw.ABDm)),std(S_raw.sf*S_raw.ABDi)./sqrt(length(S_raw.ABDi))],...
%     '-ko','LineWidth',2,'CapSize',0,'MarkerFaceColor','w','MarkerSize',9);

%%

range = [0,100];

S_2 = load('../matFiles/slopes_CO_top500_unscaled_Pot_2.mat');
S_2.Int= removeOutliers(S_2.Int_pot2,range);
S_2.Vest = removeOutliers(S_2.Vest_pot2,range);
S_2.ABDm = removeOutliers(S_2.ABDm_pot2,range);
S_2.ABDi = removeOutliers(S_2.ABDi_pot2,range);

S_2.SF = GF_mean(1)/ nanmean(S_2.Int);


S_5 = load('../matFiles/slopes_CO_top500_unscaled_Pot_5.mat');
S_5.Int= removeOutliers(S_5.Int_pot5,range);
S_5.Vest = removeOutliers(S_5.Vest_pot5,range);
S_5.ABDm = removeOutliers(S_5.ABDm_pot5,range);
S_5.ABDi = removeOutliers(S_5.ABDi_pot5,range);

S_5.SF = GF_mean(1)/ nanmean(S_5.Int);


S_10 = load('../matFiles/slopes_CO_top500_unscaled_Pot_10.mat');
S_10.Int= removeOutliers(S_10.Int_pot10,range);
S_10.Vest = removeOutliers(S_10.Vest_pot10,range);
S_10.ABDm = removeOutliers(S_10.ABDm_pot10,range);
S_10.ABDi = removeOutliers(S_10.ABDi_pot10,range);

S_10.SF = GF_mean(1)/ nanmean(S_10.Int);

cols = cbrewer('qual','Dark2',5);


subplot(2,4,1)

distributionPlot([S_raw_actual.sf*rawData(:,1);S_raw_actual.sf*rawData(:,2)],'xValues',1,'histOri','left','color','k','widthDiv',[2,1],'showMM',1,'histOpt',0,'divFactor',[-1:1:20]);
hold on
distributionPlot(S_raw_actual.sf*rawData(:,3),'xValues',2,'histOri','left','color','k','widthDiv',[2,1],'showMM',1,'histOpt',0,'divFactor',[-1:1:20]);
distributionPlot([S_raw_actual.sf*rawData(:,4);S_raw_actual.sf*rawData(:,5)],'xValues',3,'histOri','left','color','k','widthDiv',[2,1],'showMM',1,'histOpt',0,'divFactor',[-1:1:20]);
%distributionPlot(S_raw.sf*rawData(:,5),'xValues',4,'histOri','left','color','k','widthDiv',[2,1],'showMM',1,'histOpt',0,'divFactor',[-1:1:20]);

distributionPlot([model_SF*Int_rm(:)],'xValues',1,'histOri','right','color','r','widthDiv',[2,2],'showMM',1,'histOpt',0,'divFactor',[min(model_SF*Int_rm(:)):1:max(model_SF*Int_rm(:))]);
distributionPlot(model_SF*Vest_rm(:),'xValues',2,'histOri','right','color','r','widthDiv',[2,2],'showMM',1,'histOpt',0,'divFactor',[min(model_SF*Vest_rm(:)):1:max(model_SF*Vest_rm(:))]);
distributionPlot([model_SF*ABDm_rm(:); model_SF*ABDi_rm(:)],'xValues',3,'histOri','right','color','r','widthDiv',[2,2],'showMM',1,'histOpt',0,'divFactor',[min(model_SF*ABDm_rm(:)):1:max(model_SF*ABDi_rm(:))]);
%distributionPlot(S_2.SF*S_2.ABDi','xValues',4,'histOri','right','color','r','widthDiv',[2,2],'showMM',1,'histOpt',0,'divFactor',[-1:1:max(S_2.ABDi)]);
title('red - synapse assignment jitter');
set(gca,'YGrid','on','YMinorGrid','on','XTick',[1,2,3],'XTickLabels',{'VPNI','DO','ABD'},'YLim',[-5,25]);
%daspect([1,6,1]);

subplot(2,4,2)

distributionPlot([S_raw_actual.sf*rawData(:,1);S_raw_actual.sf*rawData(:,2)],'xValues',1,'histOri','left','color','k','widthDiv',[2,1],'showMM',1,'histOpt',0,'divFactor',[-1:1:20]);
hold on
distributionPlot(S_raw_actual.sf*rawData(:,3),'xValues',2,'histOri','left','color','k','widthDiv',[2,1],'showMM',1,'histOpt',0,'divFactor',[-1:1:20]);
distributionPlot([S_raw_actual.sf*rawData(:,4);S_raw_actual.sf*rawData(:,5)],'xValues',3,'histOri','left','color','k','widthDiv',[2,1],'showMM',1,'histOpt',0,'divFactor',[-1:1:20]);
%distributionPlot(S_raw.sf*rawData(:,5),'xValues',4,'histOri','left','color','k','widthDiv',[2,1],'showMM',1,'histOpt',0,'divFactor',[-1:1:20]);

distributionPlot(S_2.SF*S_2.Int','xValues',1,'histOri','right','color',cols(1,:),'widthDiv',[2,2],'showMM',1,'histOpt',0,'divFactor',[-1:1:max(S_2.Int)]);
distributionPlot(S_2.SF*S_2.Vest','xValues',2,'histOri','right','color',cols(1,:),'widthDiv',[2,2],'showMM',1,'histOpt',0,'divFactor',[-1:1:max(S_2.Vest)]);
distributionPlot([S_2.SF*S_2.ABDm';S_2.SF*S_2.ABDi'],'xValues',3,'histOri','right','color',cols(1,:),'widthDiv',[2,2],'showMM',1,'histOpt',0,'divFactor',[-1:1:max(S_2.ABDm)]);
%distributionPlot(S_2.SF*S_2.ABDi','xValues',4,'histOri','right','color','r','widthDiv',[2,2],'showMM',1,'histOpt',0,'divFactor',[-1:1:max(S_2.ABDi)]);
set(gca,'YGrid','on','YMinorGrid','on','XTick',[1,2,3],'XTickLabels',{'VPNI','DO','ABD'});
%daspect([1,6,1]);
title('red - 2\mum jitter');


subplot(2,4,3)

% distributionPlot([S_raw_actual.sf*rawData(:,1);S_raw_actual.sf*rawData(:,2)],'xValues',1,'histOri','left','color','k','widthDiv',[2,1],'showMM',1,'histOpt',0,'divFactor',[-1:1:20]);
% hold on
% distributionPlot(S_raw_actual.sf*rawData(:,3),'xValues',2,'histOri','left','color','k','widthDiv',[2,1],'showMM',1,'histOpt',0,'divFactor',[-1:1:20]);
% distributionPlot([S_raw_actual.sf*rawData(:,4);S_raw_actual.sf*rawData(:,5)],'xValues',3,'histOri','left','color','k','widthDiv',[2,1],'showMM',1,'histOpt',0,'divFactor',[-1:1:20]);
% %distributionPlot(S_raw.sf*rawData(:,5),'xValues',4,'histOri','left','color','k','widthDiv',[2,1],'showMM',1,'histOpt',0,'divFactor',[-1:1:20]);

distributionPlot([expt_SF_1_99*exptData_1_99(:,1);expt_SF_1_99*exptData_1_99(:,2)],'xValues',1,'histOri','right','color','r','widthDiv',[2,2],'showMM',1,'histOpt',0,'divFactor',[0:1:50]);
hold on;
distributionPlot(expt_SF_1_99*exptData_1_99(:,3),'xValues',2,'histOri','right','color','r','widthDiv',[2,2],'showMM',1,'histOpt',0,'divFactor',[0:1:50]);
distributionPlot([expt_SF_1_99*exptData_1_99(:,4);expt_SF_1_99*exptData_1_99(:,5)],'xValues',3,'histOri','right','color','r','widthDiv',[2,2],'showMM',1,'histOpt',0,'divFactor',[0:1:50]);


distributionPlot(S_5.SF*S_5.Int','xValues',1,'histOri','left','color',cols(1,:),'widthDiv',[2,1],'showMM',1,'histOpt',0,'divFactor',[-1:1:max(S_5.Int)]);
distributionPlot(S_5.SF*S_5.Vest','xValues',2,'histOri','left','color',cols(1,:),'widthDiv',[2,1],'showMM',1,'histOpt',0,'divFactor',[-1:1:max(S_5.Vest)]);
distributionPlot([S_5.SF*S_5.ABDm';S_5.SF*S_5.ABDi'],'xValues',3,'histOri','left','color',cols(1,:),'widthDiv',[2,1],'showMM',1,'histOpt',0,'divFactor',[-1:1:max(S_5.ABDm)]);
%distributionPlot(S_5.SF*S_5.ABDi','xValues',4,'histOri','right','color','r','widthDiv',[2,2],'showMM',1,'histOpt',0,'divFactor',[-1:1:max(S_5.ABDi)]);
title('red - 5\mum jitter');
set(gca,'YGrid','on','YMinorGrid','on','XTick',[1,2,3],'XTickLabels',{'VPNI','DO','ABD'},'YLim',[0,15]);
%daspect([1,6,1]);

subplot(2,4,4)

distributionPlot([S_raw_actual.sf*rawData(:,1);S_raw_actual.sf*rawData(:,2)],'xValues',1,'histOri','left','color','k','widthDiv',[2,1],'showMM',1,'histOpt',0,'divFactor',[-1:1:20]);
hold on
distributionPlot(S_raw_actual.sf*rawData(:,3),'xValues',2,'histOri','left','color','k','widthDiv',[2,1],'showMM',1,'histOpt',0,'divFactor',[-1:1:20]);
distributionPlot([S_raw_actual.sf*rawData(:,4);S_raw_actual.sf*rawData(:,5)],'xValues',3,'histOri','left','color','k','widthDiv',[2,1],'showMM',1,'histOpt',0,'divFactor',[-1:1:20]);
%distributionPlot(S_raw.sf*rawData(:,5),'xValues',4,'histOri','left','color','k','widthDiv',[2,1],'showMM',1,'histOpt',0,'divFactor',[-1:1:20]);

distributionPlot(S_10.SF*S_10.Int','xValues',1,'histOri','right','color',cols(3,:),'widthDiv',[2,2],'showMM',1,'histOpt',0,'divFactor',[-1:1:max(S_10.Int)]);
distributionPlot(S_10.SF*S_10.Vest','xValues',2,'histOri','right','color',cols(3,:),'widthDiv',[2,2],'showMM',1,'histOpt',0,'divFactor',[-1:1:max(S_10.Vest)]);
distributionPlot([S_10.SF*S_10.ABDm';S_10.SF*S_10.ABDi'],'xValues',3,'histOri','right','color',cols(3,:),'widthDiv',[2,2],'showMM',1,'histOpt',0,'divFactor',[-1:1:max(S_10.ABDm)]);
%distributionPlot(S_10.SF*S_10.ABDi','xValues',4,'histOri','right','color','r','widthDiv',[2,2],'showMM',1,'histOpt',0,'divFactor',[-1:1:max(S_10.ABDi)]);
title('red - 10\mum jitter');
set(gca,'YGrid','on','YMinorGrid','on','XTick',[1,2,3],'XTickLabels',{'VPNI','DO','ABD'});
%daspect([1,6,1]);


figure;

subplot(2,2,1)

errorbar([1,2,3,4],[nanmean([S_raw_actual.sf*rawData(:,1);S_raw_actual.sf*rawData(:,2)]),...
    nanmean(S_raw_actual.sf*rawData(:,3)),nanmean(S_raw_actual.sf*rawData(:,4)),nanmean(S_raw_actual.sf*rawData(:,5))],...
    [nanstd([S_raw_actual.sf*rawData(:,1);S_raw_actual.sf*rawData(:,2)]),...
    nanstd(S_raw_actual.sf*rawData(:,3)),nanstd(S_raw_actual.sf*rawData(:,4)),nanstd(S_raw_actual.sf*rawData(:,5))]./...
    [sqrt(sum(~isnan([S_raw_actual.sf*rawData(:,1);S_raw_actual.sf*rawData(:,2)]))),sqrt(sum(~isnan(S_raw_actual.sf*rawData(:,3)))),sqrt(sum(~isnan(S_raw_actual.sf*rawData(:,4)))),sqrt(sum(~isnan(S_raw_actual.sf*rawData(:,5))))],...
    '-ko','lineWidth',2);
hold on;

errorbar([1,2,3.5],[nanmean([[expt_SF_1_99*exptData_1_99(:,1);expt_SF_1_99*exptData_1_99(:,2)]]),...
    nanmean(expt_SF_1_99*exptData_1_99(:,3)),nanmean([expt_SF_1_99*exptData_1_99(:,4);expt_SF_1_99*exptData_1_99(:,5)])],...
    [nanstd([[expt_SF_1_99*exptData_1_99(:,1);expt_SF_1_99*exptData_1_99(:,2)]]),...
    nanstd(expt_SF_1_99*exptData_1_99(:,3)),nanstd([expt_SF_1_99*exptData_1_99(:,4);expt_SF_1_99*exptData_1_99(:,5)])]./...
    [sqrt(sum(~isnan([expt_SF_1_99*exptData_1_99(:,1);expt_SF_1_99*exptData_1_99(:,2)]))),...
    sqrt(sum(~isnan(exptData_1_99(:,3)))),sqrt(sum(~isnan(exptData_1_99(:,4))) + sum(~isnan(exptData_1_99(:,5))))],...
    '-o','lineWidth',2);

errorbar([1,2,3,4],[nanmean(S_2.SF*S_2.Int),nanmean(S_2.SF*S_2.Vest),...
    nanmean(S_2.SF*S_2.ABDm),nanmean(S_2.SF*S_2.ABDi)],...
    [nanstd(S_2.SF*S_2.Int),nanstd(S_2.SF*S_2.Vest),...
    nanstd(S_2.SF*S_2.ABDm),nanstd(S_2.SF*S_2.ABDi)]./...
[sqrt(sum(~isnan(S_2.Int))),sqrt(sum(~isnan(S_2.Vest))),sqrt(sum(~isnan(S_2.ABDm))),sqrt(sum(~isnan(S_2.ABDi)))],...
    '-o','lineWidth',2);

errorbar([1,2,3,4],[nanmean(S_5.SF*S_5.Int),nanmean(S_5.SF*S_5.Vest),...
    nanmean(S_5.SF*S_5.ABDm),nanmean(S_5.SF*S_5.ABDi)],...
    [nanstd(S_5.SF*S_5.Int),nanstd(S_5.SF*S_5.Vest),...
    nanstd(S_5.SF*S_5.ABDm),nanstd(S_5.SF*S_5.ABDi)]./...
    [sqrt(sum(~isnan(S_5.Int))),sqrt(sum(~isnan(S_5.Vest))),sqrt(sum(~isnan(S_5.ABDm))),sqrt(sum(~isnan(S_5.ABDi)))],...
    '-o','lineWidth',2);

errorbar([1,2,3,4],[nanmean(S_10.SF*S_10.Int),nanmean(S_10.SF*S_10.Vest),...
    nanmean(S_10.SF*S_10.ABDm),nanmean(S_10.SF*S_10.ABDi)],...
    [nanstd(S_10.SF*S_10.Int),nanstd(S_10.SF*S_10.Vest),...
    nanstd(S_10.SF*S_10.ABDm),nanstd(S_10.SF*S_10.ABDi)]./...
    [sqrt(sum(~isnan(S_10.Int))),sqrt(sum(~isnan(S_10.Vest))),sqrt(sum(~isnan(S_10.ABDm))),sqrt(sum(~isnan(S_10.ABDi)))],...
    '-o','lineWidth',2);
legend({'W_{ij}','Functional imaging','W_{ij}-2\mum','W_{ij}-5\mum','W_{ij}-10\mum'});
set(gca,'XTick',[1,2,3,4],'XTickLabels',{'VPNI','Vest','ABDm','ABDi'});
box off;
axis square;


subplot(2,2,2) % medians

errorbar([1,2,3,4],[nanmedian([S_raw_actual.sf*rawData(:,1);S_raw_actual.sf*rawData(:,2)]),...
    nanmedian(S_raw_actual.sf*rawData(:,3)),nanmedian(S_raw_actual.sf*rawData(:,4)),nanmedian(S_raw_actual.sf*rawData(:,5))],...
    [nanstd([S_raw_actual.sf*rawData(:,1);S_raw_actual.sf*rawData(:,2)]),...
    nanstd(S_raw_actual.sf*rawData(:,3)),nanstd(S_raw_actual.sf*rawData(:,4)),nanstd(S_raw_actual.sf*rawData(:,5))],...
    '-ko','lineWidth',2);
hold on;

errorbar([1,2,3.5],[nanmedian([[expt_SF_1_99*exptData_1_99(:,1);expt_SF_1_99*exptData_1_99(:,2)]]),...
    nanmedian(expt_SF_1_99*exptData_1_99(:,3)),nanmedian([expt_SF_1_99*exptData_1_99(:,4);expt_SF_1_99*exptData_1_99(:,5)])],...
    [nanstd([[expt_SF_1_99*exptData_1_99(:,1);expt_SF_1_99*exptData_1_99(:,2)]]),...
    nanstd(expt_SF_1_99*exptData_1_99(:,3)),nanstd([expt_SF_1_99*exptData_1_99(:,4);expt_SF_1_99*exptData_1_99(:,5)])],...
    '-o','lineWidth',2);

errorbar([1,2,3,4],[nanmedian(S_2.SF*S_2.Int_pot2),nanmedian(S_2.SF*S_2.Vest_pot2),...
    nanmedian(S_2.SF*S_2.ABDm_pot2),nanmedian(S_2.SF*S_2.ABDi_pot2)],...
    [nanstd(S_2.SF*S_2.Int),nanstd(S_2.SF*S_2.Vest),...
    nanstd(S_2.SF*S_2.ABDm),nanstd(S_2.SF*S_2.ABDi)],...
    '-o','lineWidth',2);

errorbar([1,2,3,4],[nanmedian(S_5.SF*S_5.Int_pot5),nanmedian(S_5.SF*S_5.Vest_pot5),...
    nanmedian(S_5.SF*S_5.ABDm_pot5),nanmedian(S_5.SF*S_5.ABDi_pot5)],...
    [nanstd(S_5.SF*S_5.Int),nanstd(S_5.SF*S_5.Vest),...
    nanstd(S_5.SF*S_5.ABDm),nanstd(S_5.SF*S_5.ABDi)],...
    '-o','lineWidth',2);

errorbar([1,2,3,4],[nanmedian(S_10.SF*S_10.Int_pot10),nanmedian(S_10.SF*S_10.Vest_pot10),...
    nanmedian(S_10.SF*S_10.ABDm_pot10),nanmedian(S_10.SF*S_10.ABDi_pot10)],...
    [nanstd(S_10.SF*S_10.Int),nanstd(S_10.SF*S_10.Vest),...
    nanstd(S_10.SF*S_10.ABDm),nanstd(S_10.SF*S_10.ABDi)],...
    '-o','lineWidth',2);
legend({'Actual','Functional','2um','5um','10um'});
set(gca,'XTick',[1,2,3,4],'XTickLabels',{'VPNI','Vest','ABDm','ABDi'});
box off;
axis square;


subplot(2,2,3)

errorbar([1,2,3],[nanmean([S_raw_actual.sf*rawData(:,1);S_raw_actual.sf*rawData(:,2)]),...
    nanmean(S_raw_actual.sf*rawData(:,3)),nanmean([S_raw_actual.sf*rawData(:,4);S_raw_actual.sf*rawData(:,5)])],...
    [nanstd([S_raw_actual.sf*rawData(:,1);S_raw_actual.sf*rawData(:,2)]),...
    nanstd(S_raw_actual.sf*rawData(:,3)),nanstd([S_raw_actual.sf*rawData(:,4);S_raw_actual.sf*rawData(:,5)])]./...
    [sqrt(sum(~isnan([S_raw_actual.sf*rawData(:,1);S_raw_actual.sf*rawData(:,2)]))),sqrt(sum(~isnan(S_raw_actual.sf*rawData(:,3)))),sqrt(sum(~isnan([S_raw_actual.sf*rawData(:,4);S_raw_actual.sf*rawData(:,5)])))],...
    '-ko','lineWidth',2);
hold on;

errorbar([1,2,3],[nanmean([[expt_SF_1_99*exptData_1_99(:,1);expt_SF_1_99*exptData_1_99(:,2)]]),...
    nanmean(expt_SF_1_99*exptData_1_99(:,3)),nanmean([expt_SF_1_99*exptData_1_99(:,4);expt_SF_1_99*exptData_1_99(:,5)])],...
    [nanstd([[expt_SF_1_99*exptData_1_99(:,1);expt_SF_1_99*exptData_1_99(:,2)]]),...
    nanstd(expt_SF_1_99*exptData_1_99(:,3)),nanstd([expt_SF_1_99*exptData_1_99(:,4);expt_SF_1_99*exptData_1_99(:,5)])]./...
    [sqrt(sum(~isnan([expt_SF_1_99*exptData_1_99(:,1);expt_SF_1_99*exptData_1_99(:,2)]))),...
    sqrt(sum(~isnan(exptData_1_99(:,3)))),sqrt(sum(~isnan(exptData_1_99(:,4))) + sum(~isnan(exptData_1_99(:,5))))],...
    '-o','lineWidth',2);

errorbar([1,2,3],[nanmean(S_2.SF*S_2.Int),nanmean(S_2.SF*S_2.Vest),...
    nanmean([S_2.SF*S_2.ABDm,S_2.SF*S_2.ABDi])],...
    [nanstd(S_2.SF*S_2.Int),nanstd(S_2.SF*S_2.Vest),...
    nanstd([S_2.SF*S_2.ABDm,S_2.SF*S_2.ABDi])]./...
[sqrt(sum(~isnan(S_2.Int))),sqrt(sum(~isnan(S_2.Vest))),sqrt(sum(~isnan([S_2.ABDm,S_2.ABDi])))],...
    '-o','lineWidth',2);

errorbar([1,2,3],[nanmean(S_5.SF*S_5.Int),nanmean(S_5.SF*S_5.Vest),...
    nanmean([S_5.SF*S_5.ABDm,S_5.SF*S_5.ABDi])],...
    [nanstd(S_5.SF*S_5.Int),nanstd(S_5.SF*S_5.Vest),...
    nanstd([S_5.SF*S_5.ABDm,S_5.SF*S_5.ABDi])]./...
    [sqrt(sum(~isnan(S_5.Int))),sqrt(sum(~isnan(S_5.Vest))),sqrt(sum(~isnan([S_5.ABDm,S_5.ABDi])))],...
    '-o','lineWidth',2);

errorbar([1,2,3],[nanmean(S_10.SF*S_10.Int),nanmean(S_10.SF*S_10.Vest),...
    nanmean([S_10.SF*S_10.ABDm,S_10.SF*S_10.ABDi])],...
    [nanstd(S_10.SF*S_10.Int),nanstd(S_10.SF*S_10.Vest),...
    nanstd([S_10.SF*S_10.ABDm,S_10.SF*S_10.ABDi])]./...
    [sqrt(sum(~isnan(S_10.Int))),sqrt(sum(~isnan(S_10.Vest))),sqrt(sum(~isnan([S_10.ABDm,S_10.ABDi])))],...
    '-o','lineWidth',2);
legend({'W_{ij}','Functional imaging','W_{ij}-2\mum','W_{ij}-5\mum','W_{ij}-10\mum'});
set(gca,'XTick',[1,2,3,4],'XTickLabels',{'VPNI','Vest','ABDm','ABDi'});
box off;
axis square;
%% box plots

figure;
subplot(1,2,1)
cols = cbrewer('qual','Dark2',5);

boxplot( S_raw_actual.sf*rawData(:),...
    [ones(numel(rawData(:,1:2)),1);2*ones(numel(rawData(:,3)),1);3*ones(numel(rawData(:,4:5)),1)],...
    'Positions',[1,2,3],'colors','k','Jitter',0.1,...
    'PlotStyle','compact');
hold on;

plot([1,2,3],[nanmean([S_raw_actual.sf*rawData(:,1);S_raw_actual.sf*rawData(:,2)]),...
    nanmean(S_raw_actual.sf*rawData(:,3)),nanmean([S_raw_actual.sf*rawData(:,4);S_raw_actual.sf*rawData(:,5)])],...
    'w^','MarkerSize',7)

boxplot(expt_SF_1_99*exptData_1_99(:),...
    [1*ones(numel(exptData_1_99(:,1:2)),1);2*ones(numel(exptData_1_99(:,3)),1);3*ones(numel(exptData_1_99(:,4:5)),1)],...
    'Positions',[1.1,2.1,3.1],'colors','r','Jitter',0.1,...
    'PlotStyle','compact');

plot([1.1,2.1,3.1],[nanmean([expt_SF_1_99*exptData_1_99(:,1);expt_SF_1_99*exptData_1_99(:,2)]),...
    nanmean(expt_SF_1_99*exptData_1_99(:,3)),nanmean([expt_SF_1_99*exptData_1_99(:,4);expt_SF_1_99*exptData_1_99(:,5)])],...
    'w^','MarkerSize',7);

boxplot([S_2.SF*S_2.Int' ;S_2.SF*S_2.Vest'; S_2.SF*S_2.ABDm' ;S_2.SF*S_2.ABDi'],...
    [1*ones(numel(S_2.Int),1);2*ones(numel(S_2.Vest),1);3*[ones(numel(S_2.ABDm),1);ones(numel(S_2.ABDi),1)]],...
    'Positions',[1.2,2.2,3.2],'colors',cols(1,:),'Jitter',0.1,...
    'PlotStyle','compact');

plot([1.2,2.2,3.2],[nanmean(S_2.SF*S_2.Int),nanmean(S_2.SF*S_2.Vest),...
    nanmean([S_2.SF*S_2.ABDm,S_2.SF*S_2.ABDi])],...
    'w^','MarkerSize',7);

boxplot([S_5.SF*S_5.Int' ;S_5.SF*S_5.Vest'; S_5.SF*S_5.ABDm' ;S_5.SF*S_5.ABDi'],...
    [1*ones(numel(S_5.Int),1);2*ones(numel(S_5.Vest),1);3*[ones(numel(S_5.ABDm),1);ones(numel(S_5.ABDi),1)]],...
    'Positions',[1.3,2.3,3.3],'colors',cols(2,:),'Jitter',0.1,...
    'PlotStyle','compact');

plot([1.3,2.3,3.3],[nanmean(S_5.SF*S_5.Int),nanmean(S_5.SF*S_5.Vest),...
    nanmean([S_5.SF*S_5.ABDm,S_5.SF*S_5.ABDi])],...
    'w^','MarkerSize',7);


boxplot([S_10.SF*S_10.Int' ;S_10.SF*S_10.Vest'; S_10.SF*S_10.ABDm' ;S_10.SF*S_10.ABDi'],...
    [1*ones(numel(S_10.Int),1);2*ones(numel(S_10.Vest),1);3*[ones(numel(S_10.ABDm),1);ones(numel(S_10.ABDi),1)]],...
    'Positions',[1.4,2.4,3.4],'colors',cols(3,:),'Jitter',0.1,...
    'PlotStyle','compact');

plot([1.4,2.4,3.4],[nanmean(S_10.SF*S_10.Int),nanmean(S_10.SF*S_10.Vest),...
    nanmean([S_10.SF*S_10.ABDm,S_10.SF*S_10.ABDi])],...
     'w^','MarkerSize',7);

set(gca,'YLim',[-2,16],'YGrid','on','YMinorGrid','on');
box off;

%%

figure;
subplot(2,2,4)

distributionPlot([S_raw_actual.sf*rawData(:,1);S_raw_actual.sf*rawData(:,2)],'xValues',1,'histOri','left','color','k','widthDiv',[2,1],'showMM',1,'histOpt',0,'divFactor',[0:1:50]);
hold on
distributionPlot(S_raw_actual.sf*rawData(:,3),'xValues',2,'histOri','left','color','k','widthDiv',[2,1],'showMM',1,'histOpt',0,'divFactor',[0:1:50]);
distributionPlot([S_raw_actual.sf*rawData(:,4);S_raw_actual.sf*rawData(:,5)],'xValues',3,'histOri','left','color','k','widthDiv',[2,1],'showMM',1,'histOpt',0,'divFactor',[0:1:50]);

distributionPlot([expt_SF_1_99*exptData_1_99(:,1);expt_SF_1_99*exptData_1_99(:,2)],'xValues',1,'histOri','right','color','r','widthDiv',[2,2],'showMM',1,'histOpt',0,'divFactor',[0:1:50]);
distributionPlot(expt_SF_1_99*exptData_1_99(:,3),'xValues',2,'histOri','right','color','r','widthDiv',[2,2],'showMM',1,'histOpt',0,'divFactor',[0:1:50]);
distributionPlot([expt_SF_1_99*exptData_1_99(:,4);expt_SF_1_99*exptData_1_99(:,5)],'xValues',3,'histOri','right','color','r','widthDiv',[2,2],'showMM',1,'histOpt',0,'divFactor',[0:1:50]);

set(gca,'YLim',[-5,20],'YGrid','on','YMinorGrid','on','XTick',[1,2,3],'XTickLabels',{'VPNI','DO','ABD'});
daspect([1,6,1]);


figure;


subplot(2,2,1)
errorbar([1.1,2,3,4],macaque_SF*macaque_mean,macaque_SF*macaque_SE,'-o','LineWidth',2,'CapSize',0,'MarkerFaceColor','w','MarkerSize',9);
hold on;
errorbar([0.9,2,3,4.1],cat_SF*cat_mean,cat_SF*cat_SE,'-mo','LineWidth',2,'CapSize',0,'MarkerFaceColor','w','MarkerSize',9);
errorbar([1,2,3,4],GF_mean,GF_SE,'-o','LineWidth',2,'CapSize',0,'MarkerFaceColor','w','MarkerSize',9);
errorbar([mean([S_raw_actual.sf*S_raw_actual.r46,S_raw_actual.sf*S_raw_actual.r78]),mean(S_raw_actual.sf*S_raw_actual.Vest),mean(S_raw_actual.sf*S_raw_actual.ABDm),mean(S_raw_actual.sf*S_raw_actual.ABDi)],...
    [std(S_raw_actual.sf*S_raw_actual.Int)./sqrt(length(S_raw_actual.Int)),std(S_raw_actual.sf*S_raw_actual.Vest)./sqrt(length(S_raw_actual.Vest)),std(S_raw_actual.sf*S_raw_actual.ABDm)./sqrt(length(S_raw_actual.ABDm)),std(S_raw_actual.sf*S_raw_actual.ABDi)./sqrt(length(S_raw_actual.ABDi))],...
    '-ko','LineWidth',2,'CapSize',0,'MarkerFaceColor','w','MarkerSize',9);


errorbar([1.05,2,3,4],[nanmean([model_SF*modelData_1000(:,1);model_SF*modelData_1000(:,2)]),...
    nanmean(model_SF*modelData_1000(:,3)),...
    nanmean(model_SF*modelData_1000(:,4)),...
    nanmean(model_SF*modelData_1000(:,5))],...
    [nanstd([model_SF*modelData_1000(:,1);model_SF*modelData_1000(:,2)])./sqrt(sum(~isnan([model_SF*modelData_1000(:,1);model_SF*modelData_1000(:,2)]))),...
    nanstd(model_SF*modelData_1000(:,3))./sqrt(sum(~isnan(modelData_1000(:,3)))),...
    nanstd(model_SF*modelData_1000(:,4))./sqrt(sum(~isnan(modelData_1000(:,4)))),...
    nanstd(model_SF*modelData_1000(:,5))./sqrt(sum(~isnan(modelData_1000(:,5))))],...
    '-o','lineWidth',2,'CapSize',0,'MarkerFaceColor','w','MarkerSize',9);


errorbar([0.95,2,3.5],[nanmean([expt_SF_1_99*exptData_1_99(:,1);expt_SF_1_99*exptData_1_99(:,2)]),...
    nanmean(expt_SF_1_99*exptData_1_99(:,3)),...
    nanmean([expt_SF_1_99*exptData_1_99(:,4);expt_SF_1_99*exptData_1_99(:,5)])],...
    [nanstd([expt_SF_1_99*exptData_1_99(:,1);expt_SF_1_99*exptData_1_99(:,2)])./sqrt(sum(~isnan([expt_SF_1_99*exptData_1_99(:,1);expt_SF_1_99*exptData_1_99(:,2)]))),...
    nanstd(expt_SF_1_99*exptData_1_99(:,3))./sqrt(sum(~isnan(exptData_1_99(:,3)))),...
    nanstd([expt_SF_1_99*exptData_1_99(:,4);expt_SF_1_99*exptData_1_99(:,5)])./...
    (sqrt(sum(~isnan(exptData_1_99(:,4))) + sum(~isnan(exptData_1_99(:,5)))))],...
    '-or','lineWidth',2,'CapSize',0,'MarkerFaceColor','w','MarkerSize',9);

box off;
axis square;
set(gca,'XTick',[1,2,3,4],'XTickLabels',{'VPNI','DO','ABDm','ABDi'});
ylabel('K')
offsetAxes(gca);
legend({'mac','cat','gf','model_{raw}','model_{1000}','expt'},'Location','bestoutside');

%% Kvalues by top 50 bins

fileName = 'ranking_noMirror_cleaned_squared.csv';
fid = fopen(fileName);
% cellID,fanOut,fanIn,toIntegrator,fromInteg,toMotor,class,centrality
SpectralWithRS.complete = textscan(fid,'%d %d %d %d %d %d %d %s %f','Delimiter',',');
fclose(fid);


load MatOrder_CO_top500_2blocks_gamma038_08062020.mat
load cellIDType_CO_top500_2blocks_gamma038_08062020.mat
cellIDType = cellstr(cellIDType_CO_top500_2blocks_gamma038_08062020);

[loA,loB] = ismember(MatOrder_CO_top500_2blocks_gamma038_08062020,SpectralWithRS.complete{2});
Int_locs = find(strcmp(cellIDType,'_Int_'));
Int_ranks = loB(strcmp(cellIDType,'_Int_'));

%top50 bins

for i = 1:11
    temp  = Int_locs(Int_ranks<50*i);
    noInts(i) = size(temp,1);
    A(i,:) = [temp(end),max(Int_ranks)];
end

%%


function[population_rm] = removeOutliers(population,range)
% ptiles = prctile(population,[range(1),range(2)],1);
% low = population>population(ptiles(1,:));
% high =  population<ptiles(2,:);
% temp = population.*low;
% population_rm = temp.*high;

%population = population(:);
low = prctile(population,range(1));
high = prctile(population,range(2));
if size(population,2) == 1
    population_rm = population((population>low) & (population<high));
else
    index = (population>low) & (population<high);
    population_rm = population.*index;
    population_rm(population_rm==0) = NaN;
    
end
end

