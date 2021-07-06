function renderLandmarksOnImageGrid(true_landmarks, estimated_landmarks, ncols, one_image_size, ROI_centroids, ROI_half_size, col_true, col_estimated, landmarked_bones, MarkerSize, LineWidth, DS)
num_frames = size(estimated_landmarks{1}{1},1);
num_views = size(estimated_landmarks,1);
for view=1:num_views
    for frame=1:num_frames
        total_indx = (frame-1)*num_views + view-1; %(k-1)*num_views+(view-1);
        img_grid_indx_y = floor(total_indx/ncols);
        img_grid_indx_x = total_indx-img_grid_indx_y*ncols;
        landmark_rendering_offset = [img_grid_indx_x img_grid_indx_y].*one_image_size;
        for k=1:length(landmarked_bones)
            if(~isempty(ROI_centroids))
                ROI_offset = ROI_centroids((frame-1)*num_views+view,:)-ROI_half_size;
            else
                ROI_offset = [0 0];
            end
            for i=1:size(estimated_landmarks{view,landmarked_bones(k)},1)
                if(size(estimated_landmarks{view,landmarked_bones(k)}{i},3)==0), continue; end
                for j=1:size(estimated_landmarks{view,landmarked_bones(k)}{i},3)
                    plot(estimated_landmarks{view,landmarked_bones(k)}{i}(frame,1,j)/DS-ROI_offset(1)+landmark_rendering_offset(1), estimated_landmarks{view,landmarked_bones(k)}{i}(frame,2,j)/DS-ROI_offset(2)+landmark_rendering_offset(2), '+', 'Color', col_estimated(i,:), 'MarkerSize', MarkerSize, 'LineWidth', LineWidth);
                    if(~isempty(true_landmarks))
                        plot(true_landmarks{view,landmarked_bones(k)}.pos2D{i}(frame,1,j)/DS-ROI_offset(1)+landmark_rendering_offset(1), true_landmarks{view,landmarked_bones(k)}.pos2D{i}(frame,2,j)/DS-ROI_offset(2)+landmark_rendering_offset(2), '+', 'Color', col_true(i,:), 'MarkerSize', MarkerSize, 'LineWidth', LineWidth);
                    end
                end
            end
        end
    end
end
