function cell_neighbors = SomasegIndex(cellOrigins1,cellOrigins2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ncells1 = size(cellOrigins1,1);
ncells2 = size(cellOrigins2,1);

isCell1 = [ones(ncells1,1);zeros(ncells2,1)];
isCell2 = [zeros(ncells1,1);ones(ncells2,1)];

allSoma = [cellOrigins1;cellOrigins2];


for i = 1:size(allSoma,1)
    dist = pdist2(allSoma(i,:),allSoma);
    [a,b] = sort(dist);
    m(i) = a(2);
    n(i) = b(2);
    clear dist;
    clear a;
    clear b;   
end

cells1_11 = sum(n(1:ncells1)<=ncells1);
cells1_12 = sum(n(1:ncells1)>ncells1);

cells2_21 = sum(n(ncells1+1:end)<=ncells1);
cells2_22 = sum(n(ncells1+1:end)>ncells1);

cell_neighbors = [cells1_11,cells1_12;cells2_21,cells2_22];

end

