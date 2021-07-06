function noisy_label = removeLabelNoise(noisy_label, threshold)
num_label = max(noisy_label(:));
background_label = 0;
if(~exist('threshold','var')), threshold = 0.05; end

noise_voxels = cell(num_label,1);
for i=1:num_label, noise_voxels{i} = false(size(noisy_label)); end
for i=1:num_label
    fprintf('removing noise of label %d/%d...', i, num_label); timerID = tic(); 
    props = regionprops3(noisy_label==i, 'Volume', 'VoxelIdxList');
    noise_islands = props.Volume < sum(props.Volume)*threshold;
    noisy_label(cell2mat(props.VoxelIdxList(noise_islands))) = background_label;
    fprintf('done in %.3f sec\n', toc(timerID));
%     [L, NUM] = bwlabeln(noisy_label==i);
%     island_volume = zeros(NUM,1);
%     for j=1:NUM
%         island_volume(j) = sum(L(:)==j);
%     end
%     island_volume_ratio = island_volume / sum(island_volume);
%     for j=1:NUM
%         if(island_volume_ratio(j)<threshold)
%             noise_voxels{i}(L==j) = true;
%         end
%     end
end
% for i=1:num_label
%     noisy_label(noise_voxels{i}) = background_label;
% end
