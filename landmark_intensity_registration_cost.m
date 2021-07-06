function [cost, trans_param_4x4_image_coord, trans_param_4x4_local_coord] = landmark_intensity_registration_cost(param, opt_data, visualization_on, force_visualize)
    if(~exist('force_visualize','var')), force_visualize = false; end
    P = size(param,2);    % Population size

    if(isempty(opt_data.landmark_tbls))
        num_frames = 1;
    else
%         num_transforms = size(opt_data.landmark_tbls,2);
        num_frames = size(opt_data.landmark_tbls{1}.pos2D{1},1);
    end
    if(opt_data.CPU_par)
        num_parallel = min(opt_data.num_transforms*6*P, maxNumCompThreads);
    else
        num_parallel = 1;
    end
    
    [trans_param_4x4_image_coord, ~, ~, trans_param_4x4_local_coord] = ...
        ComputeTransformationParameters(param, opt_data.num_ctrl_pnts, opt_data.num_transforms, num_frames, num_parallel, opt_data.ctrl_scale, ...
        opt_data.offset_scale, opt_data.global_transform_6xMxN, opt_data.volume_offset_4x4xN(:,:,opt_data.target_object), opt_data.fix_ctrl_offset, opt_data.rotation_center_local_4x4xN(:,:,opt_data.target_object));

    % landmark cost
    if(isempty(opt_data.landmark_tbls) || opt_data.cost_alpha==1.0)
        cost_landmark = zeros(1,P);
    else
        mean_reprojection_error = ProjectLandmarks(opt_data.landmark_tbls, trans_param_4x4_image_coord, opt_data.ProjectionMatrices_pix, false);
        cost_landmark = opt_data.cost_landmark_scaling * max(squeeze(mean(mean(mean_reprojection_error,1,'omitnan'),2,'omitnan'))'-opt_data.landmark_cost_threshold, 0);
    end
    
    % rigidity regularization (minimize norm of transformation relative to
    % the reference bone. the first index bone is set to the reference bone
    % for now)
    if(opt_data.num_transforms>1 && ~isempty(opt_data.rigidity_constraint_object) &&  opt_data.lambda_rigidity_regularization>0)
        rigidity_regularization = opt_data.lambda_rigidity_regularization * compute_relative_transformation_norm(trans_param_4x4_image_coord(:,:,opt_data.rigidity_constraint_object,:,:)); %reshape(mean(mean(var(trans_param_NxF_array_reshaped(:,opt_data.rigidity_constraint_object,:,:),0,2),4),1),1,P);
    else
        rigidity_regularization = zeros(1,P);
    end
    
    % smoothness regularization
    if(opt_data.lambda_smoothness_regularization==0 || num_frames<=2)
        smoothness_regularization = zeros(1,P);
    else
        smoothness_regularization = opt_data.lambda_smoothness_regularization * compute_transformation_curvature(trans_param_4x4_image_coord); % reshape(mean(mean(mean(abs(curvature),4),2),1),1,P);
    end
    
    % joint regularization
    if(~isempty(opt_data.joint_constraint_tbl) && opt_data.lambda_joint_regularization>0)
        all_joint_regularization = zeros(height(opt_data.joint_constraint_tbl),P);
        for i=1:height(opt_data.joint_constraint_tbl)
            switch opt_data.joint_constraint_tbl.constraint_type{i}
                case 'spherical'
                    trans_param_4x4_image_coord_4Nx4_base = reshape(permute(trans_param_4x4_local_coord(:,:,opt_data.joint_constraint_tbl.base_object(i),:,:), [1 3 4 5 2]),[],4);
                    joint_center_base_global = reshape(trans_param_4x4_image_coord_4Nx4_base*[opt_data.joint_constraint_tbl.coordinate_4x4_wrt_base{i}(1:3,4); 1], 4, P, num_frames);
                    trans_param_4x4_image_coord_4Nx4_ref = reshape(permute(trans_param_4x4_local_coord(:,:,opt_data.joint_constraint_tbl.ref_object(i),:,:), [1 3 4 5 2]),[],4);
                    joint_center_ref_global = reshape(trans_param_4x4_image_coord_4Nx4_ref*[opt_data.joint_constraint_tbl.coordinate_4x4_wrt_ref{i}(1:3,4); 1], 4, P, num_frames);
                    euclidean_dist = sqrt(sum((joint_center_base_global(1:3,:,:)-joint_center_ref_global(1:3,:,:)).^2,1));
                    all_joint_regularization(i,:) = mean(max(euclidean_dist-opt_data.joint_constraint_tbl.constraint_threshold(i),0),3);
            end
        end
        joint_regularization = opt_data.lambda_joint_regularization * mean(all_joint_regularization,1);
    else
        joint_regularization = zeros(1,P);
    end
    
    % image similarity
    numGPUs = length(opt_data.GPU_IDs);
    if(opt_data.enable_image && opt_data.cost_alpha>0)
        num_views = size(opt_data.ProjectionMatrices_pix,1)/3;
        if(P<numGPUs)  % for multi-GPU environment (RegTools complains error if the number of image set is not a multiple of the number of GPUs)
            trans_param_4x4_local_coord = repmat(trans_param_4x4_local_coord,[1 1 1 numGPUs 1]);
            GetMergedFloatingImages( opt_data.regTools, opt_data.volumePlans, permute(trans_param_4x4_local_coord(:,:,:,:,opt_data.image_similarity_frame),[1 2 3 5 4]), num_views, opt_data.target_object);
%         fprintf('rendering finished: %.2f MB available on the GPU\n', opt_data.regTools.GPUmemCheck/1024/1024);
            cost_similarity = -opt_data.regTools.ComputeSimilarityMeasure( opt_data.similarity_measure_plan_id, opt_data.SimilarityMeasureType, size(trans_param_4x4_local_coord,4) )'; % cost: row vector
            cost_similarity = cost_similarity(:,1:P);
        else
            indx = array_split_1D(1:P, max(1,floor(P/opt_data.maxParallelRendering)));
            cost_similarity = zeros(1,size(trans_param_4x4_local_coord,4));
            for i=1:length(indx)
                if(isempty(indx{i})), continue; end
                GetMergedFloatingImages( opt_data.regTools, opt_data.volumePlans, permute(trans_param_4x4_local_coord(:,:,:,indx{i},opt_data.image_similarity_frame),[1 2 3 5 4]), num_views, opt_data.target_object);
                cost_similarity(1,indx{i}) = -opt_data.regTools.ComputeSimilarityMeasure( opt_data.similarity_measure_plan_id, opt_data.SimilarityMeasureType, length(indx{i}) )'; % cost: row vector
            end
        end
        cost_similarity = cost_similarity*opt_data.cost_similarity_scaling;
        if(0)   % for debugging
            debug_rendering_cost;
        end
    else
        cost_similarity = zeros(1,P);
    end
    
    cost = (1-opt_data.cost_alpha)*cost_landmark + opt_data.cost_alpha*cost_similarity + rigidity_regularization + smoothness_regularization + joint_regularization;
    opt_data.iteration_count = opt_data.iteration_count + P;
    opt_data.cost_log = [opt_data.cost_log; [cost' cost_landmark' cost_similarity' rigidity_regularization' smoothness_regularization' joint_regularization']];
    elapsed_time = toc(opt_data.TimerID);
    
    % output current cost for monitoring
    [~, min_cost_indx] = min(cost, [], 2);
    fprintf('%s, %d iter, landmark:%.3f pix, similarity:%.3f, rigidity:%.3f, smoothness:%.3f, joint:%.3f, cost:%f, %.3f sec, %.3f ite/sec (%d ctrl pts, %d DoF, %d param, PopSize: %d)\n', ...
           opt_data.experiment_ID, opt_data.iteration_count, cost_landmark(min_cost_indx), cost_similarity(min_cost_indx), rigidity_regularization(min_cost_indx), smoothness_regularization(min_cost_indx), joint_regularization(min_cost_indx), cost(min_cost_indx), ...
           elapsed_time, (opt_data.iteration_count-opt_data.previous_iteration)/(elapsed_time-opt_data.previous_elapsed_time), ...
           opt_data.num_ctrl_pnts, opt_data.num_transforms*6, size(param,1), opt_data.PopSize);
       
    opt_data.previous_elapsed_time = elapsed_time;
    opt_data.previous_iteration = opt_data.iteration_count;

    if(visualization_on)
        if(opt_data.iteration_count == 1 || opt_data.iteration_count>(opt_data.previous_show_iteration + opt_data.progress_show_step) || force_visualize)
            [min_cost, min_cost_indx] = min(cost, [], 2);
            min_cost_indx = min_cost_indx(1);   % for a special case where multiple (exactly the same) minimums were found
            [trans_param_4x4_image_coord, trans_param_NxF, scaled_param, trans_param_4x4_local_coord] = ComputeTransformationParameters(param(:,min_cost_indx), opt_data.num_ctrl_pnts, opt_data.num_transforms, num_frames, 1, ...
                opt_data.ctrl_scale, opt_data.offset_scale, opt_data.global_transform_6xMxN, opt_data.volume_offset_4x4xN(:,:,opt_data.target_object), opt_data.fix_ctrl_offset, opt_data.rotation_center_local_4x4xN(:,:,opt_data.target_object));
            if(~isempty(opt_data.landmark_tbls))
                [mean_reprojection_error, all_projected] = ProjectLandmarks(opt_data.landmark_tbls, trans_param_4x4_image_coord, opt_data.ProjectionMatrices_pix, true);         
            else
                mean_reprojection_error = []; all_projected = [];
            end
%             trans_param_NxF_array_reshaped = reshape(RegTools.convert4x4ToTransRot_multi(reshape(trans_param_4x4_image_coord,4,4,[])),6,opt_data.num_transforms, 1, num_frames);

            show_progress(opt_data, trans_param_4x4_local_coord, trans_param_NxF, scaled_param, all_projected, mean_reprojection_error, min_cost, size(param,1));
        end
    end
end
