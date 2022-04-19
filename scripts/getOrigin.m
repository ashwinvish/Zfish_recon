function [origin] = getOrigin(cellID)
% *getOrigin* returns a 3x1 =, x,y,z coordinated of cell in cellID
% in Z-brian space
% cellID is the Id of the cell

for i = 1:numel(cellID)
    if (isExistReRoot(cellID(i)) == 1)
        tree = SwctoZbrian(cellID(i));
        origin(i,:) = [tree{1}.X(1), tree{1}.Y(1), tree{1}.Z(1)];
        clear tree;
    else
        origin(i,:) = [0,0,0];
    end
end
end