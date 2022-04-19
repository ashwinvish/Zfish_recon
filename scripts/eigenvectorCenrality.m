% eigenvector centrality

% clc;
% clear;

load ConnMatrixPre_cleaned.mat
load AllCells_cleaned.mat
load cellIDs.mat;


format longG

bilateral = 0;
cellClass = cell(size(AllCells_cleaned));
names = fieldnames(cellIDs);

for i = 1:size(names,1)
    tempIDs = getfield(cellIDs,char(names(i)));
    locs = ismember(AllCells_cleaned,tempIDs);
    cellClass(locs) = names(i);
    clear tempIDs;
    clear locs;
end
cellClass(cellfun(@isempty,cellClass)) = {'NaN'};

W = ConnMatrixPre_cleaned;

% manual cleanup
toClean = [77150,77658 81172 82194]; % ABD neurons
%%%% if RS neurons are zeroed
%toClean = [toClean,76202,77099,76562,77259,77260,77265,77267,77268,77441,77449,77456,77694,77695,77931,78566,78577,78940,79244,79395,79961,80327,81410,81611,82217,82218,82220]; % RS neurons
%%%%

toCleanLocs = find(ismember(AllCells_cleaned,toClean));
W(:,toCleanLocs) = 0;

% IDs of r78contra dendrites (unique). Each dendrite should have no outgoing synapses
dendrites = [76199, 76200, 76182, 76183, 76185, 76186, 76189, 76191, 76188, 77582, 77605, 79040, 76399, 76828, 76829, 76826, 76289, 76542, 76832, 76838, 76877];
% best matches with orphan axons (some are repeated). Each axon should have no incoming synapses
axons = [78687, 78651, 76666, 80219, 76666, 76666, 78677, 79950, 77869, 80219, 79134, 78923, 76675, 78903, 81682, 80241, 80248, 78615, 77773, 78615, 80242];

%  synthesize bilateral model
%  mirror so there are twice as many segments
%  glue together axons and dendrites from opposite sides (correspondences above)
%  orphan axons should have zero incoming synapses to start with
%  this will remove their outgoing synapses too


Wipsi = W;
Crossing = zeros(size(W)); % synthesize connections crossing the midline
for i = 1:length(axons)
    dendriteind = find(AllCells_cleaned == dendrites(i));
    axonind = find(AllCells_cleaned == axons(i));
    %transfer outgoing connections from axon to corresponding dendrite on other side
    Crossing(:, dendriteind) = Wipsi(:, axonind);
    %zero out outgoing connections from axon
    Wipsi(:, axonind) = 0;
end

%This is the 2n x 2n connection matrix for the bilateral model
Wbilateral = [Wipsi Crossing; Crossing Wipsi];

%principal eigenvector is first one since eigenvalues are sorted

if bilateral ==1
    
    [vl,~] = eigs(Wbilateral',1,'largestreal');
    [vr,~] = eigs(Wbilateral,1,'largestreal');
    
    vl = -vl;
    vr = -vr; % to match the julia convention
    
    
    centrality_squared = sqrt(abs(vl).*abs(vr));
    n = size(W, 1) ; % number of neurons without mirroring

    [~,sortInd_squared] = sort(centrality_squared(1:n),'descend');
    
    centrality_prod = vl.*vr;
    
    outDegree = sum(Wbilateral,1);
    inDegree = sum(Wbilateral,2);
    
    centrality_deg = sqrt(outDegree.*inDegree');
    
    for i = 1:length(centrality_prod)
        if centrality_prod(i) < -1e-6
            centrality_prod(i) = -centrality_prod(i);
        end
    end
    
    [~,sortInd_prod] = sort(centrality_prod(1:n),'descend');
    [~,sortInd_deg] = sort(centrality_deg,'descend');

        
    filename_squared = 'ranking_Mirror_cleaned_noRS_squared.csv';
    filename_prod = 'ranking_Mirror_cleaned_noRS_prod.csv';
    filename_deg = 'ranking_Mirror_cleaned_noRS_deg.csv';

else
    
    [vl,~] = eigs(W',1,'largestreal');
    [vr,~] = eigs(W,1,'largestreal');
    
    vl = -vl;
    vr = -vr; % to match the julia convention
    
    
    centrality_squared = sqrt(abs(vl).*abs(vr));
    [~,sortInd_squared] = sort(centrality_squared,'descend');
    
    centrality_prod = vl.*vr;
    
    outDegree = sum(W,1);
    inDegree = sum(W,2);
    
    centrality_deg = sqrt(outDegree.*inDegree');
    
    
    for i = 1:length(centrality_prod)
        if centrality_prod(i) < -1e-6
            centrality_prod(i) = -centrality_prod(i);
        end
    end
    
    n = size(W, 1) ; % number of neurons without mirroring
    [~,sortInd_prod] = sort(centrality_prod,'descend');
    
    [~,sortInd_deg] = sort(centrality_deg,'descend');
    
     filename_squared = 'ranking_noMirror_cleaned_squared.csv';
     filename_prod = 'ranking_noMirror_cleaned_prod.csv';
    filename_deg = 'ranking_noMirror_cleaned_deg.csv';

end
   

%% plots


subplot(4,4,1)
plot(abs(vl(sortInd_squared(1:500))));
hold on
plot(abs(vr(sortInd_squared(1:500))));
plot(centrality_squared(sortInd_squared(1:500)),'-k','LineWidth',2),
set(gca,'YScale','log')
box off;
axis square;
legend({'left centrality','right centrality','mean'})

subplot(4,4,2)
histogram(W(W>0),'FaceColor','none','EdgeColor','k','LineWidth',2);
set(gca,'YScale','log');
axis square;
box off;
%offsetAxes(gca);

subplot(4,4,3)
plot(outDegree(sortInd_squared));
set(gca,'YScale','log');

hold on
plot(inDegree(sortInd_squared));
set(gca,'YScale','log');

box off;

EC_cells = AllCells_cleaned(sortInd_squared);
DC_cells = AllCells_cleaned(sortInd_deg);
[~,temp2] = ismember(EC_cells,DC_cells);

subplot(4,4,4)
scatter(1:540,temp2(1:540),3,'ko');
box off;
f = showfit(ezfit(1:540,temp2(1:540),'linear'),'fitcolor','r','dispeqboxmode','off');
axis square;
set(gca,'XLim',[0,600],'YLim',[0,600],'XTick',[0,200,400,600]);
text(650,650,sprintf('r = %0.3f',f.r));
xlabel('Eigen centrality rank');
ylabel('Degree centrality rank')


%%
    
    fileID = fopen(filename_squared,'w');
    for i = 1:n
        fprintf(fileID, '%d,%d,%d,%d,%d,%d,%d,%s,%e\n',sortInd_squared(i),AllCells_cleaned(sortInd_squared(i)),sum(W(:,sortInd_squared(i)),1)',...
            sum(W(sortInd_squared(i),:),2),0,0,0,string(cellClass(sortInd_squared(i))),centrality_squared(sortInd_squared(i)));
    end
    fclose(fileID);
    
    
    fileID = fopen(filename_prod,'w');
    for i = 1:n
        fprintf(fileID, '%d,%d,%d,%d,%d,%d,%d,%s,%e\n',sortInd_prod(i),AllCells_cleaned(sortInd_prod(i)),sum(W(:,sortInd_prod(i)),1)',...
            sum(W(sortInd_prod(i),:),2),0,0,0,string(cellClass(sortInd_prod(i))),centrality_prod(sortInd_prod(i)));
    end
    fclose(fileID);

    
    fileID = fopen(filename_deg,'w');
    for i = 1:n
        fprintf(fileID, '%d,%d,%d,%d,%d,%d,%d,%s,%e\n',sortInd_deg(i),AllCells_cleaned(sortInd_deg(i)),sum(W(:,sortInd_deg(i)),1)',...
            sum(W(sortInd_deg(i),:),2),0,0,0,string(cellClass(sortInd_deg(i))),centrality_deg(sortInd_deg(i)));
    end
    fclose(fileID);   

