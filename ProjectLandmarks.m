function [mean_reprojection_error, all_projected] = ProjectLandmarks(landmark_tbls, trans_param_4x4, ProjectionMatrices_pix, return_projected_pnts)
    % project landmarks and calculate mean reprojection error
    [~, ~, num_transforms, P, num_frames] = size(trans_param_4x4);
    num_views = size(ProjectionMatrices_pix,1)/3;
    mean_reprojection_error = zeros(num_views, num_transforms, P);
    all_PM = ProjectionMatrices_pix * reshape(trans_param_4x4,4,[]);
    all_PM_realigned = reshape(permute(reshape(all_PM,3,num_views,4,num_transforms,P,[]),[1 6 5 3 4 2]),[],4,num_transforms,num_views);
    if(return_projected_pnts)
        all_projected = cell(num_views, num_transforms);
    end
    for j=1:num_views           % loop over views
        for k=1:num_transforms     % loop over number of transformations
            pos3D_linear = cell2mat(landmark_tbls{j,k}.pos3D);
            p = all_PM_realigned(:,:,k,j) * [pos3D_linear'; ones(1,size(pos3D_linear,1))];
            if(return_projected_pnts)
                projected_xy_con = reshape( p([1:3:end 2:3:end],:)./repmat(p(3:3:end,:),2,1), num_frames, 2, []);
                projected_xy = permute(projected_xy_con, [1 3 2] );
                nl = cell2mat(cellfun(@size, landmark_tbls{j,k}.pos3D, 'UniformOutput', false));
                all_projected{j,k} = squeeze(mat2cell(projected_xy_con, num_frames, 2, nl(:,1)));
            else
                projected_xy = permute( reshape( p([1:3:end 2:3:end],:)./repmat(p(3:3:end,:),2,1), num_frames*P, 2, []), [1 3 2] );
            end
            landmark_2D_linear = repmat(permute(cell2mat(cellfun(@(x) permute(x,[3 1 2]),landmark_tbls{j,k}.pos2D,'UniformOutput',false)),[2 1 3]),P,1,1);
            mean_reprojection_error(j,k,:) = mean(mean(reshape(sqrt(sum( (landmark_2D_linear(1:num_frames*P,:,:)-projected_xy).^2, 3)),num_frames,P,[]),1),3);
        end
    end

