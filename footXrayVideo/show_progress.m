function show_progress(opt_data, trans_param_4x4_local_coord, trans_param_NxF, scaled_param, all_projected, mean_reprojection_error, min_cost,  num_parameters)
    [hs,vs,tb,bb,lb,rb] = deal(0.07, 0.07, 0.12, 0.1, 0.03, 0.03);
    clf;
    colormap(gray(256));
    colorbar_title = {'Fixed landmarks', 'Floating landmarks'};
    num_bones_each_transform = zeros(opt_data.num_transforms,1);
    all_bone_names = {};
    fixed_clim = opt_data.fixed_clim; grad_clim = opt_data.grad_clim; moving_clim = opt_data.moving_clim;
    for j=1:opt_data.num_transforms
        num_bones_each_transform(j) = length(opt_data.landmark_tbls{1,j}.bone_name);
        all_bone_names = cat(1, all_bone_names, opt_data.landmark_tbls{1,j}.bone_name);
    end
    if(opt_data.num_transforms==1)
        landmarked_bones = find(cell2mat(cellfun(@(x) ~isempty(x), opt_data.landmark_tbls{1}.pos2D, 'UniformOutput', false)));
    else
        landmarked_bones = find(cell2mat(cellfun(@(x) ~isempty(x.pos2D{1}), opt_data.landmark_tbls(1,:), 'UniformOutput', false)));
    end
    num_bones = sum(num_bones_each_transform);
    col_true = hsv(length(landmarked_bones));
    col_estimated = min(col_true+0.6,1); %parula(num_bones);
    colorbar_col = cat(3, col_true, col_estimated);
    hs2 = 0.005; vs2 = 0.02; tb2 = 0.07; bb2 = 0.03; lb2 = 0.01;
    if(length(opt_data.registration_frame)>1)
        ncols = 2;
    else
        ncols = 1;
    end
    if(opt_data.enable_image)
        num_views = size(opt_data.ProjectionMatrices_pix,1)/3;
        trans_param_min_cost = repmat(trans_param_4x4_local_coord(:,:,:,1,opt_data.image_similarity_frame), [1 1 1 1 length(opt_data.GPU_IDs)]);
        floating = GetMergedFloatingImages( opt_data.regTools, opt_data.volumePlans, trans_param_min_cost, num_views, opt_data.target_object );
        floating = floating(:,:,1:length(opt_data.image_similarity_frame)*num_views);
        cost_similarity = opt_data.regTools.ComputeSimilarityMeasure( opt_data.similarity_measure_plan_id, opt_data.SimilarityMeasureType, length(opt_data.GPU_IDs) );
        fixed = opt_data.regTools.GetSimilarityMeasureComputationPlan(opt_data.similarity_measure_plan_id, 0, [], [], 1);
        intermediate_image = opt_data.regTools.GetSimilarityMeasureComputationPlan(opt_data.similarity_measure_plan_id, 1, [], [], 1);
        landmark_centroids_view1 = floor(opt_data.landmark_centroids_2D{1}(opt_data.image_similarity_frame,:));
        landmark_centroids_view2 = floor(opt_data.landmark_centroids_2D{2}(opt_data.image_similarity_frame,:));
        landmark_centroids = reshape([[landmark_centroids_view1(:,1)'; landmark_centroids_view2(:,1)'] [landmark_centroids_view1(:,2)'; landmark_centroids_view2(:,2)']],[],2);
        ROI_size = opt_data.zoom_ROI_size;
        ROI_fixed = zeros(ROI_size*2+1, ROI_size*2+1, size(fixed,3));
        ROI_floating = zeros(ROI_size*2+1, ROI_size*2+1, size(floating,3));
        ROI_intermediate_image = zeros(ROI_size*2+1, ROI_size*2+1, size(intermediate_image,3));
        ROI_centroids = zeros(size(landmark_centroids));
        for i=1:size(fixed,3)
            ROI_centroids(i,:) = [max(landmark_centroids(i,1),1+ROI_size) max(landmark_centroids(i,2),1+ROI_size)]; 
            ROI_centroids(i,:) = [min(ROI_centroids(i,1),size(fixed,1)-ROI_size) min(ROI_centroids(i,2),size(fixed,2)-ROI_size)]; 
            ROI_fixed(:,:,i) = fixed(ROI_centroids(i,1)+(-ROI_size:ROI_size), ROI_centroids(i,2)+(-ROI_size:ROI_size),i);
            fixed(ROI_centroids(i,1) + [-ROI_size+[0 1 2] +ROI_size+[-2 -1 0]], (ROI_centroids(i,2)-ROI_size):(ROI_centroids(i,2)+ROI_size), i) = max(fixed(:));
            fixed((ROI_centroids(i,1)-ROI_size):(ROI_centroids(i,1)+ROI_size), ROI_centroids(i,2) + [-ROI_size+[0 1 2] +ROI_size+[-2 -1 0]], i) = max(fixed(:));
            ROI_floating(:,:,i) = floating(ROI_centroids(i,1)+(-ROI_size:ROI_size), ROI_centroids(i,2)+(-ROI_size:ROI_size),i);
            floating(ROI_centroids(i,1) + [-ROI_size+[0 1 2] +ROI_size+[-2 -1 0]], (ROI_centroids(i,2)-ROI_size):(ROI_centroids(i,2)+ROI_size), i) = max(floating(:));
            floating((ROI_centroids(i,1)-ROI_size):(ROI_centroids(i,1)+ROI_size), ROI_centroids(i,2) + [-ROI_size+[0 1 2] +ROI_size+[-2 -1 0]], i) = max(floating(:));
            ROI_intermediate_image(:,:,i) = intermediate_image(ROI_centroids(i,1)+(-ROI_size:ROI_size), ROI_centroids(i,2)+(-ROI_size:ROI_size),i);
            intermediate_image(ROI_centroids(i,1) + [-ROI_size+[0 1 2] +ROI_size+[-2 -1 0]], (ROI_centroids(i,2)-ROI_size):(ROI_centroids(i,2)+ROI_size), i) = max(intermediate_image(:));
            intermediate_image((ROI_centroids(i,1)-ROI_size):(ROI_centroids(i,1)+ROI_size), ROI_centroids(i,2) + [-ROI_size+[0 1 2] +ROI_size+[-2 -1 0]], i) = max(intermediate_image(:));
        end
        if(opt_data.LCN_sigma>0)
            fixed_clim = [-2 2];    % for LCN
        end
        ROI_fixed_rgb = permute( repmat( mat2gray(ROI_fixed,fixed_clim), [1 1 1 3]), [1 2 4 3]);
        for i=1:size(ROI_fixed_rgb,4)
            contour = 255 - canny('image',uint8(ROI_floating(:,:,i)./max(ROI_floating(:)).*255), 'thigh', 0.5);
            contour = imdilate(contour, strel('square', 2));
            ROI_fixed_rgb(:,:,1,i) = ROI_fixed_rgb(:,:,1,i) + double(logical(contour));
        end

        label_slice = reshape([cellstr(strcat(num2str(opt_data.registration_frame'),'-l')) cellstr(strcat(num2str(opt_data.registration_frame'),'-o'))]',[],1);
        msubplot(1,7,1,hs2,vs2,tb2,bb2,lb2,rb); hold on;
        im(permute(fixed,[2 1 3]),'ncols',ncols,'labelslice',label_slice); set(gca,'clim',fixed_clim);  title('Fixed images', 'FontSize', 12); axis off;
        renderLandmarksOnImageGrid(opt_data.landmark_tbls, all_projected, ncols, [size(fixed,1) size(fixed,2)], [], [], col_true, col_estimated, landmarked_bones, 5, 1, opt_data.DS_2D); set(gca,'ydir','reverse');
        msubplot(1,7,2,hs2,vs2,tb2,bb2,lb2,rb); hold on;
        imc(permute(ROI_fixed_rgb,[2 1 3 4]),'ncols',ncols,'labelslice',label_slice); set(gca,'clim',fixed_clim);   title('ROI', 'FontSize', 12); axis off;
        renderLandmarksOnImageGrid(opt_data.landmark_tbls, all_projected, ncols, [1 1]*ROI_size*2+3, ROI_centroids, ROI_size, col_true, col_estimated, landmarked_bones, 5, 1, opt_data.DS_2D); set(gca,'ydir','reverse');
        msubplot(1,7,3,hs2,vs2,tb2,bb2,lb2,rb); hold on;
        im(permute(floating,[2 1 3]),'ncols',ncols,'labelslice',label_slice); set(gca,'clim',moving_clim); title('Floating images', 'FontSize', 12); axis off;
        renderLandmarksOnImageGrid([], all_projected, ncols, [size(fixed,1) size(fixed,2)], [], [], col_true, col_estimated, landmarked_bones, 5, 1, opt_data.DS_2D); set(gca,'ydir','reverse');
        msubplot(1,7,4,hs2,vs2,tb2,bb2,lb2,rb); hold on;
        im(permute(ROI_floating,[2 1 3]),'ncols',ncols,'labelslice',label_slice); set(gca,'clim',moving_clim); title('ROI', 'FontSize', 12); axis off;
        renderLandmarksOnImageGrid([], all_projected, ncols, [1 1]*ROI_size*2+3, ROI_centroids, ROI_size, col_true, col_estimated, landmarked_bones, 5, 1, opt_data.DS_2D); set(gca,'ydir','reverse');
        msubplot(1,7,5,hs2,vs2,tb2,bb2,lb2,rb); im(permute(intermediate_image,[2 1 3]),'ncols',ncols,'labelslice',label_slice); set(gca,'clim',grad_clim/1000); title('Gradient Correlation', 'FontSize', 12); axis off;
        msubplot(1,7,6,hs2,vs2,tb2,bb2,lb2,rb); im(permute(ROI_intermediate_image,[2 1 3]),'ncols',ncols,'labelslice',label_slice); set(gca,'clim',grad_clim/1000); title('ROI', 'FontSize', 12); axis off;
        nrows = 3; indx_offset = 24; 
        colorbar_pos = [0.14 0.06 0.2 0.01; 0.35 0.06 0.2 0.01];
%                 for i=1:2
%                     msubplot(1,9,6,hs2,vs2,tb2,bb2,lb2,rb);
%                     colormap(gca,colorbar_col(:,:,i));
% %                     set(gca,'clim',[0 num_bones]);
%                     c = colorbar('Location','South','Position',colorbar_pos(i,:),'TickLabels',all_bone_names,'Ticks',0.5:1:num_bones-0.5);
%                     c.Label.String = colorbar_title{i};
%                     c.Label.FontSize = 10;
%                 end
    else
        nrows = 1; indx_offset = 0; 
        colorbar_pos = [0.14 0.06 0.2 0.01; 0.35 0.06 0.2 0.01];
    end
    for i=1:2
        if(opt_data.enable_image), continue; end
        if(opt_data.enable_image)
            msubplot(nrows,9,indx_offset+i,0.03,0.03,tb,0.05,lb,rb);
        else
            msubplot(nrows,3,indx_offset+i,0.03,0.03,tb,0.05,lb,rb);
        end
        hold on;
        for k=1:opt_data.num_transforms
            if(opt_data.num_transforms==1)
                renderLandmarkTrajectory(opt_data.landmark_tbls{i,k}.pos2D, all_projected{i,k}, opt_data.DimSize_2D(i,:), col_true, col_estimated, landmarked_bones);
            else
                renderLandmarkTrajectory(opt_data.landmark_tbls{i,k}.pos2D, all_projected{i,k}, opt_data.DimSize_2D(i,:), col_true(k,:), col_estimated(k,:), landmarked_bones);
            end
        end
        title({ sprintf('%s view',opt_data.view_direction{i}), 'landmark trajectory', ...
                'reprojection error:', sprintf('%.3f pix',mean(mean_reprojection_error(i,:),'omitnan'))}, 'FontSize', 12);
        colormap(gca,colorbar_col(:,:,i));
        set(gca,'clim',[0 length(landmarked_bones)]);
        c = colorbar('Location','South','Position',colorbar_pos(i,:),'TickLabels',all_bone_names(landmarked_bones),'Ticks',0.5:1:num_bones-0.5);
        c.Label.String = colorbar_title{i};
        c.Label.FontSize = 10;
    end

%             offset = opt_data.global_transform_6xF + p((1:num_transforms*6) + num_transforms*opt_data.num_ctrl_pnts*6,min_cost_indx).*opt_data.offset_scale;
    if(opt_data.num_transforms==1)
        msubplot(2,3,3,hs,vs,tb,bb,lb,rb);
        plot_trans_param(trans_param_NxF, scaled_param, opt_data.global_transform_6xMxN(:,:,1));
        title('Transformation parameters (zero mean)');
%     elseif(num_transforms==4)
%         indx = [5 6 11 12];
%         for i=1:4
%             msubplot(3,6,indx(i),0.05,vs,tb,bb,lb,rb);
%             plot_trans_param(trans_param_NxF((1:6)+(i-1)*6,:), scaled_param((1:6)+(i-1)*6,:), opt_data.global_transform_6xMxN(:,:,i));
%             title(sprintf('%s transform (zero mean)',opt_data.landmark_tbls{1,i}.bone_name{1}));
%         end
    end

    if(length(opt_data.cost_log)>1)
        if(opt_data.num_transforms==1)
            msubplot(2,3,6,hs,vs,tb,bb,lb,rb);
            FontSize = 10;
        else
%                     msubplot(3,9,27,0.03,0.03,tb,0.05,lb,0.02);
            msubplot(1,7,7,hs2+0.02,vs2+0.02,tb2+0.02,bb2+0.02,lb2+0.05,rb); 
            FontSize = 8;
        end
        plot(1:size(opt_data.cost_log,1),opt_data.cost_log)
        set(gca, 'FontSize', 10, 'xlim', [1 length(opt_data.cost_log)]); %, 'ylim', [0 50]);
%                 ylabel('mean reprojection error [mm]'); 
        legend('cost', 'landmark error', 'image similarity', sprintf('rigidity (\\lambda: %.1e)', opt_data.lambda_rigidity_regularization), ...
            sprintf('smoothness (\\lambda: %.1e)', opt_data.lambda_smoothness_regularization), ...
            sprintf('joint (\\lambda: %.1e)', opt_data.lambda_joint_regularization), 'FontSize', FontSize);
        xlabel('# of iterations');
    end

    btitle(sprintf('%s, Iteration: %d, cost: %f, %.3f sec (%d control points, %d DoF, %d parameters, PopSize: %d, alpha: %.2f)', ...
                opt_data.experiment_ID, opt_data.iteration_count, min_cost, ...
                toc(opt_data.TimerID), opt_data.num_ctrl_pnts, opt_data.num_transforms*6, num_parameters, opt_data.PopSize, opt_data.cost_alpha), 18, [0 0 0], 'none', 0.97);
    drawnow;
    if(~isempty(opt_data.writerObj))
        writeVideo(opt_data.writerObj, getframe(gcf));
    end
    opt_data.previous_show_iteration = opt_data.iteration_count;

    if(opt_data.enable_output_debug_renderings)
        output_debug_renderings;
    end
end

function plot_trans_param(interpolated_param, scaled_param, offset)
    hold on;
    num_frames = size(interpolated_param,2);
    if(num_frames==1)
        offset_interpolated = offset;
    else
        x = linspace(0,1,size(offset,2));
        xq = linspace(0,1,num_frames);
        F = griddedInterpolant(x, 1:size(offset,2), 'spline');
        offset_interpolated = zeros(size(interpolated_param,1), num_frames);
        for i=1:size(interpolated_param,1)
            F.Values = offset(i,:);
            offset_interpolated(i,:) = F(xq);
        end
    end
    mean_trans = mean(offset_interpolated,2);
    trans_param_NxF_offset = interpolated_param; %-mean_trans*ones(1,num_frames);
    scaled_param_zeromean = scaled_param; %-mean_trans*ones(1,size(scaled_param,2));
    col = [1 0 1; 0 1 1; 1 0 0; 0 1 0; 0 0 1; 0 0 0];
%     for i=1:size(trans_param_NxF_offset,1)
%         plot(1:num_frames, trans_param_NxF_offset(i,:)', 'Color', col(i,:));
%     end
    for i=1:size(scaled_param_zeromean,1)
        plot(1:num_frames, trans_param_NxF_offset(i,:)', 'Color', col(i,:));
    end
    for i=1:size(scaled_param_zeromean,1)
        plot(linspace(1,num_frames,size(scaled_param,2)), scaled_param_zeromean(i,:)', '.', 'Color', col(i,:), 'MarkerSize', 15);
    end
    if(num_frames>1)
        set(gca,'xlim',[1 num_frames],'FontSize',9);
    end
    legend_str = {'X trans','Y trans','Z trans','X rot','Y rot','Z rot'};
    mean_str = cellfun(@(x) sprintf('(\\mu: %.2f)',x), num2cell(mean_trans),'UniformOutput',false);
    legend(strcat(legend_str', mean_str),'FontSize',6);
end

function renderLandmarkTrajectory(true_landmarks2D_cell, estimated_landmarks2D, DimSize_2D, col_true, col_estimated, landmarked_bones)
    for i=1:length(landmarked_bones) %size(true_landmarks2D_cell,1)
        for j=1:size(true_landmarks2D_cell{landmarked_bones(i)},3)
            plot(true_landmarks2D_cell{landmarked_bones(i)}(:,1,j), true_landmarks2D_cell{landmarked_bones(i)}(:,2,j), '-', 'Color', col_true(i,:));
            plot(estimated_landmarks2D{landmarked_bones(i)}(:,1,j), estimated_landmarks2D{landmarked_bones(i)}(:,2,j), '.', 'Color', col_estimated(i,:), 'MarkerSize', 5);
        end
        plot(squeeze(true_landmarks2D_cell{landmarked_bones(i)}(1,1,[1:end])), squeeze(true_landmarks2D_cell{landmarked_bones(i)}(1,2,[1:end])), '--', 'Color', col_true(i,:), 'LineWidth', 2);
        plot(squeeze(estimated_landmarks2D{landmarked_bones(i)}(1,1,[1:end])), squeeze(estimated_landmarks2D{landmarked_bones(i)}(1,2,[1:end])), '--', 'Color', col_estimated(i,:), 'LineWidth', 2);
    end
    axis equal;
    set(gca,'xlim',[0 DimSize_2D(1)],'ylim',[0 DimSize_2D(2)], 'ydir', 'reverse', 'box', 'on');
end
