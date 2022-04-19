function [ocularSelectivity,motorSum,interSum] = OSI(cellIDs,df)
% OSI ocular selectivity Index; defined as the ratio M-I/M+I
% cellIDs is the ID of the cells
% ocularSelectivity is the OSI ranges from -1 to 1
% motorSum,interSum returns the respsective sums

%format longE
motorDistribution = isMotor(cellIDs,df);
motorSum = motorDistribution(:,2)+motorDistribution(:,3);
motorSum = double(motorSum);
interSum = motorDistribution(:,4)+motorDistribution(:,5);
interSum = double(interSum);
ocularSelectivity = (motorSum-interSum)./(motorSum+interSum);
end

