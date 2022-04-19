function [edges,motorSums,interSums] = isMotor(cellID,df)
% isMotor checks to see if the cellID belong to abducens motor or
% intenuclear neurons
% edges return the number of synapses onto each neuron
% motorSums returns sum of synapses onto all motor neurons
% intersSums returns sum of synapses onto all internuclear neurons


if ~isempty(cellID)
    
    ABDr_CellIDs = [82140,82145,82143,77648,82146,77302,77301, 82194,77705,77710,77305,77672,77300,77709,77661];
    ABDc_CellIDs = [82213,77654,77646,82212,77295,81172,77657,82197,77682,77154, 77652,77658,77628,77292,77688,82195,77296];
    ABDIr_CellIDs = [78553,77150,77886,78547,77631,77158,77665,78552,77668,77618,77634,78556];
    ABDIc_CellIDs = [78574,79051,77148,77625,77144,77692,77641,79066,77640,77643];
    
    
    ABD_All = [ABDr_CellIDs,ABDc_CellIDs,ABDIr_CellIDs,ABDIc_CellIDs];
    % give each type an ID
    ABD_All_type = [ones(1,size(ABDr_CellIDs,2)),2*ones(1,size(ABDc_CellIDs,2)),...
        3*ones(1,size(ABDIr_CellIDs,2)),4*ones(1,size(ABDIc_CellIDs,2))] ;
    
    %UniqueCells = unique(cellID);
    
    for i =1:size(cellID,1)
        [A,B] = SynapticPartners(cellID(i),2,df);
        [lia, lob] =  ismember(A,ABD_All);
        temp = ABD_All_type(lob(lob>0))';
        edges(i,:) = [cellID(i),histcounts(temp,[0.9,1.9,2.9,3.1,4])];
    end
    
    motorSums =sum(edges(:,2:3),2);
    interSums = sum(edges(:,4:5),2);
    
else
    edges = [];
end

end