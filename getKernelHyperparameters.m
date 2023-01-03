function var = getKernelHyperparameters(cellIndex,treeDepth,numChild,maxGrid,minGrid,alpha)

j = treeDepth;
a = numChild;
b = a^j; % num of divisions per side

Size = maxGrid+2-minGrid;
cellSizeX1 = Size(1)/b;
LBX1 = mod(cellIndex,b)*cellSizeX1 +minGrid(1)-1;
UBX1 = (mod(cellIndex,b)+1)*cellSizeX1 +minGrid(2)-1;

cellSizeX2 = Size(2)/b;
LBX2 = mod(floor(cellIndex/b),b)*cellSizeX2 +minGrid(2)-1;
UBX2 = (mod(floor(cellIndex/b),b)+1)*cellSizeX2+minGrid(2)-1;

var = sqrt(abs(max([UBX1-LBX1,UBX2-LBX2])));


end