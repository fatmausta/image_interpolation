function [ minDist ] = RMSE3D(autoLabel, refLabel,Res)
% RMSE calculation of two 3D binary volumes using a KD-tree search
%


display('Calculating RMSE distance...');

[f1, autoSet] = isosurface(autoLabel,0);
[f2, refSet] = isosurface(refLabel,0);

tree = kdtree_build(refSet);

pntidx = kdtree_nearest_neighbor(tree, autoSet);

pntval = [refSet(pntidx,1) refSet(pntidx,2)  refSet(pntidx,3)];

minDist = sqrt((Res^2)*(autoSet(:,1)-pntval(:,1)).^2 + (Res^2)*(autoSet(:,2)-pntval(:,2)).^2 + (Res^2) *(autoSet(:,3)-pntval(:,3)).^2);


end

