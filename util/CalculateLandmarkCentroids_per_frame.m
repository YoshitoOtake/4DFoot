function landmark_centroids_2D = CalculateLandmarkCentroids_per_frame(landmark_tbls, target_object, global_transform_6xMxN, volume_offset_wrt_camera1, ProjectionMatrices_pix)

landmark_centroids_2D = cell(size(landmark_tbls,1),1);
for i=1:size(landmark_tbls,1)          % loop over views
    if(~isempty(landmark_tbls{i,1}.pos3D{1}))   % check landmark of the first bone
        all_landmarks = [];
        for j=1:length(target_object)
            all_landmarks = cat(3, all_landmarks, landmark_tbls{i,target_object(j)}.pos2D{1});
        end
        landmark_centroids_2D{i} = mean(all_landmarks,3);
    else
        num_frames = size(global_transform_6xMxN,2);
        all_centroids = zeros(num_frames,2,length(target_object));
        for j=1:num_frames
            for k=1:length(target_object)     % loop over number of bones
                projected_xyw = ProjectionMatrices_pix(:,:,i) * RegTools.convertTransRotTo4x4(global_transform_6xMxN(:,j,k)) * volume_offset_wrt_camera1(:,:,k) * [0 0 0 1]';
                all_centroids(j,:,k) = projected_xyw(1:2)./[projected_xyw(3); projected_xyw(3)];
            end
        end
        landmark_centroids_2D{i} = mean(all_centroids,3);
    end
end

