function [similarity_per_frame, all_floating_images] = calculate_similarity_per_frame(opt_data, similarity_measure_mask, img_2D_size, DownsampleRatio, img_2D_DS, param)
param = param(:,1);

GI_Sigma = 1.0;
num_frames = size(opt_data.landmark_tbls{1}.pos2D{1},1);
num_parallel = 1;
num_views = size(opt_data.ProjectionMatrices_pix,1)/3;
numGPUs = length(opt_data.GPU_IDs);
num_image_sets = length(opt_data.GPU_IDs);

[~, ~, ~, trans_param_4x4_local_coord] = ...
    ComputeTransformationParameters(param, opt_data.num_ctrl_pnts, opt_data.num_transforms, num_frames, num_parallel, opt_data.ctrl_scale, ...
    opt_data.offset_scale, opt_data.global_transform_6xMxN, opt_data.volume_offset_4x4xN(:,:,opt_data.target_object), opt_data.fix_ctrl_offset, opt_data.rotation_center_local_4x4xN(:,:,opt_data.target_object));
if(numGPUs>1)
    % for multi-GPU environment (RegTools complains error if the number of image set is not a multiple of the number of GPUs)
    trans_param_4x4_local_coord = repmat(trans_param_4x4_local_coord,[1 1 1 numGPUs 1]);
end

% calculate average similarity for confirmation purpose only
mask = opt_data.regTools.Downsample2DProjections(DownsampleRatio, reshape(similarity_measure_mask(:,:,:,opt_data.registration_frame),img_2D_size(1),img_2D_size(2),[]));
similarity_measure_plan_id = opt_data.regTools.CreateSimilarityMeasureComputationPlan( img_2D_DS, GI_Sigma, mask, num_image_sets, [], 0, 0, 0, 0, [], opt_data.SimilarityMeasureType );
GetMergedFloatingImages( opt_data.regTools, opt_data.volumePlans, permute(trans_param_4x4_local_coord(:,:,:,:,opt_data.image_similarity_frame),[1 2 3 5 4]), num_views, opt_data.target_object);
average_similarity = opt_data.regTools.ComputeSimilarityMeasure( similarity_measure_plan_id, opt_data.SimilarityMeasureType, size(trans_param_4x4_local_coord,4) ); % cost: row vector
fprintf('average similarity: %.10f\n', average_similarity(1));
opt_data.regTools.DeleteSimilarityMeasureComputationPlan( similarity_measure_plan_id );

% calculate similarity per frame
similarity_per_frame = zeros(num_frames,1,'single');
all_floating_images = zeros(img_2D_size, 'single');
fprintf('similarity per frame:');
for i=1:num_frames
    mask = opt_data.regTools.Downsample2DProjections(DownsampleRatio, reshape(similarity_measure_mask(:,:,:,opt_data.registration_frame(i)),img_2D_size(1),img_2D_size(2),[]));
    similarity_measure_plan_id = opt_data.regTools.CreateSimilarityMeasureComputationPlan( img_2D_DS(:,:,(1:2)+(opt_data.image_similarity_frame(i)-1)*2), GI_Sigma, mask, num_image_sets, [], 0, 0, 0, 0, [], opt_data.SimilarityMeasureType );
    floating_images_temp = GetMergedFloatingImages( opt_data.regTools, opt_data.volumePlans, permute(trans_param_4x4_local_coord(:,:,:,:,opt_data.image_similarity_frame(i)),[1 2 3 5 4]), num_views, opt_data.target_object);
    all_floating_images(:,:,:,i) = floating_images_temp(:,:,1:num_views);
    similarity_temp = opt_data.regTools.ComputeSimilarityMeasure( similarity_measure_plan_id, opt_data.SimilarityMeasureType, size(trans_param_4x4_local_coord,4) );
    similarity_per_frame(i) = similarity_temp(1);
    opt_data.regTools.DeleteSimilarityMeasureComputationPlan( similarity_measure_plan_id );
    fprintf(' %.10f',similarity_per_frame(i));
end
fprintf('\n');
all_floating_images = reshape(all_floating_images, img_2D_size(1), img_2D_size(2), []);

