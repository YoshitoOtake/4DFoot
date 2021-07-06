function surface_patch = prepareSurfacePatches(masked_volumes_cell, ElementSpacing_1x3)
% note: we use "smoothpatch" function here

fprintf('generating surface patches...'); timerID = tic;
num_objects = numel(masked_volumes_cell);
surface_patch = cell(num_objects, 1);
for i=1:num_objects
    fprintf('(%d/%d) done, ', i, num_objects);
    one_volume = permute(masked_volumes_cell{i},[2 1 3])>0;
    one_volume([1 end],:,:) = 0;     one_volume(:,[1 end],:) = 0;    one_volume(:,:,[1 end]) = 0;
    fv = isosurface(one_volume, 0.5);    % permute for consistency between the image coordinate and isosurfae coordinate
    fv.vertices = (fv.vertices-1-size(masked_volumes_cell{i})/2)*diag(ElementSpacing_1x3);
    surface_patch{i} = smoothpatch(reducepatch(fv,0.1));
%     surface_patch{i} = reducepatch(smoothpatch(fv,1,3), 0.1);
end
fprintf('done in %f sec\n', toc(timerID));
