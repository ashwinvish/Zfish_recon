% clc;
% clear;

load ('../matFiles/modularityScaling_500_OM_0_01_cleaned.mat');

scale = modularityScaling_500_OM_0_01_cleaned;



%% calculate variation of information

A = nchoosek(1:200,2);

for n = 1:19
    for l = 1:size(A,1)
        ARI(l,n) = rand_index(scale.clusterOrganization{n}(:,A(l,1)),scale.clusterOrganization{n}(:,A(l,2)),'adjusted');
        VarOfInf(l,n) = VOI(scale.clusterOrganization{n}(:,A(l,1)),scale.clusterOrganization{n}(:,A(l,2)));
    end
end

clear A;

scale.VOI = VarOfInf;
scale.ARI = ARI;

%% plots

%gamma = 0.1:0.05:1;
gamma = 0.1:0.05:1

subplot(2,2,1)
errorbar(gamma,cellfun(@nanmean,scale.uniqueClusters),cellfun(@nanstd,scale.uniqueClusters),'-o','LineWidth',2)
ylabel('number of modules');
title('CO-OM-top500-EigenCentrality');
box off;

subplot(2,2,3)
errorbar(gamma,mean(scale.VOI),std(scale.VOI),'-bo','LineWidth',2);
hold on
errorbar(gamma,nanmean(scale.ARI),nanstd(scale.ARI),'-ro','LineWidth',2);
errorbar(gamma,nanmean(scale.coorelations),nanstd(scale.coorelations),'-ko','LineWidth',2);
legend({'VOI','ARI','Corr'},'Location','west');
xlabel('resolution parameter \gamma')
box off;

subplot(2,2,[2,4])
violinPlot(scale.ARI,'histOri','right','color','r','widthDiv',[2,2],'showMM',1,'histOpt',1)
hold on;
violinPlot(scale.VOI,'histOri','left','color','b','widthDiv',[2 1],'showMM',1,'histOpt',1)
set(gca,'XTickLabels',gamma,'XTickLabelRotation',45);
xlabel('resolution parameter \gamma')

%%
subplot(2,3,1)
errorbar(gamma,mean(modularityScaling_600_cleaned.VOI),std(modularityScaling_600_cleaned.VOI),'-ro','LineWidth',2);
hold on
errorbar(gamma,mean(modularityScaling_500.VOI),std(modularityScaling_500.VOI),'-bo','LineWidth',2);
line([0.38,0.38],[0,1.4])
title('VOI');


subplot(2,3,2)
errorbar(gamma,mean(modularityScaling_600_cleaned.ARI),std(modularityScaling_600_cleaned.ARI),'-ro','LineWidth',2);
hold on
errorbar(gamma,mean(modularityScaling_500.ARI),std(modularityScaling_500.ARI),'-bo','LineWidth',2);
title('ARI')

subplot(2,3,3)
errorbar(gamma,mean(modularityScaling_600_cleaned.coorelations),std(modularityScaling_600_cleaned.coorelations),'-ro','LineWidth',2);
hold on
errorbar(gamma,mean(modularityScaling_500.coorelations),std(modularityScaling_500.coorelations),'-bo','LineWidth',2);
title('Corr')
legend({'DegreeCentrality','EigenCentrality'});
xlabel('resolution parameter \gamma')

