gamma = 0.1:0.05:1;
nTrials = 200;
ARI = zeros(nTrials,length(gamma));
VarInfo = zeros(nTrials,length(gamma));

for nn = 1:length(gamma)
    sprintf('gamma = %d', gamma(nn))
    cluster = [];
    meanMod = [];
    D = [];
    Dr = [];
    T =[];
    Tr = [];
    
    parfor kk = 1:nTrials
        kk;
        Q0 = -1; Q1 = 0;                            % initialize modularity values
        i =1 ;
        while abs(Q1-Q0)>1e-5                        % while modularity increases
            Q0 = Q1;                                  % perform community detection
            [M, Q1] = community_louvain(oculoMotor_syn,gamma(nn),1:size(oculoMotor_syn,1));
            i = i+1;
        end
        meanMod(kk,1) = Q1;
        cluster(:,kk) = M;
        uniqClust(1,kk) = length(unique(M));
    end
    
    cluster_reOrdered = zeros(size(cluster));
    
    for i = 1:size(cluster,2)
        unq = unique(cluster(:,i));
        locs = [];
        axialSum = [];
        motorSum = [];
        interSum = [];
        
        for j = 1:length(unq)
            [~,indMotorSum,indInterSum] = isMotor(AllCells(cluster(:,i)==unq(j)),df_cleaned);
            motorSum(j) = sum(indMotorSum);
            interSum(j) = sum(indInterSum);
            axialSum(j) = sum(isAxial(AllCells(cluster(:,i)==unq(j)),df_cleaned));
            locs{j}= find(cluster(:,i)==unq(j));
            clear indMotorSum;
            clear indInterSum;
        end
        
        % sorting order; axial to motor
        
        if length(unq) == 1
            ind = 1;
        elseif length(unq) == 2
            [~,ind] = sort(motorSum+interSum,'ascend'); % axial, motor // inter, motor
        elseif length(unq) ==3
            [~,indA] = sort(axialSum,'ascend');
            [~,indBD] = sort(motorSum+interSum,'ascend');
            ind = [indA(3),indBD(3),indBD(2)]; % axial, dual, motor // dual/none,inter, motor
        elseif length(unq) == 4
            [~,indA] = sort(axialSum,'descend');
            [~,indB] = sort(motorSum,'ascend');
            [~,indD] = sort(interSum,'ascend');
            ind = [indA(1),indD(4),indB(4),indA(2)];  % axial, inter, motor, dual/none
        elseif length(unq) == 5
            [~,indA] = sort(axialSum,'descend');
            [~,indBD] = sort(motorSum+interSum,'ascend');
            [~,indB] = sort(motorSum,'ascend');
            [~,indD] = sort(interSum,'ascend');
            
            if indBD(5) == indB(5)
                ind = [indA(1),indBD(4),indBD(5),indBD(2:3)];
            else
                ind = [indA(1),indD(5),indB(5),indBD(2:3)];
            end
        elseif length(unq) == 6
            [~,indA] = sort(axialSum,'descend');
            [~,indBD] = sort(motorSum+interSum,'ascend');
            [~,indB] = sort(motorSum,'ascend');
            [~,indD] = sort(interSum,'ascend');
            
            if indBD(6) == indB(6)
                ind = [indA(1),indBD(5),indBD(6),indBD(2:4)];
            else
                ind = [indA(1),indD(6),indB(6),indBD(2:4)];
            end
        else
            [~,indA] = sort(axialSum,'descend');
            [~,indBD] = sort(motorSum+interSum,'ascend');
            [~,indB] = sort(motorSum,'ascend');
            [~,indD] = sort(interSum,'ascend');
            ind = indB;
        end
        
        % [~,ind]  =  sort(b+d,'ascend');
        
        for j = 1:length(unq)
            cluster_reOrdered(locs{ind(j)},i) = j;
        end
        
        clear locs;
        clear ind;
        clear b;
        clear c;
        clear d;
        
    end
    
    A = nchoosek(1:nTrials,2);
    
    parfor i = 1:size(A,1)
        ARI(i,nn) = rand_index(cluster_reOrdered(:,A(i,1)),cluster_reOrdered(:,A(i,2)),'adjusted');
        VarInfo(i,nn) = VOI(cluster_reOrdered(:,A(i,1)),cluster_reOrdered(:,A(i,2)));
    end
    
    clear A;
    
    corrs = corr(cluster_reOrdered,cluster_reOrdered);
    
    corls(:,nn) = corrs(:); % all correlations
    corrsOrganization{nn} = corrs; % coorelation matrix
    clustOrganization{nn} = cluster_reOrdered; % cluster matric
    clusterOriginal{nn} = cluster;
    uniqueClusters{nn} = uniqClust; % numner of unique clusters
    modu{nn} = meanMod;
    
    clear corrs;
    
end

modularityScaling_500_OM_0_01_cleaned.coorelations = corls;
modularityScaling_500_OM_0_01_cleaned.ARI = ARI;
modularityScaling_500_OM_0_01_cleaned.VOI = VarInfo;
modularityScaling_500_OM_0_01_cleaned.coorelationMatrix = corrsOrganization;
modularityScaling_500_OM_0_01_cleaned.clusterOrganization = clustOrganization;
modularityScaling_500_OM_0_01_cleaned.clusterOriginal = clusterOriginal;
modularityScaling_500_OM_0_01_cleaned.uniqueClusters = uniqueClusters;
modularityScaling_500_OM_0_01_cleaned.modulatiry = modu;


save('modularityScaling_500_OM_0_01_cleaned.mat','modularityScaling_500_OM_0_01_cleaned');