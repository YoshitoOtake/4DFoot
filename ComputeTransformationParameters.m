function [trans_param_4x4_image_coord, trans_param_NxF_array, scaled_param, trans_param_4x4_local_coord] = ...
    ComputeTransformationParameters(p, num_ctrl_pnts, num_transforms, num_frames, num_parallel, ctrl_scale, offset_scale, global_transform_6xMxN, volume_offset_4x4xN, fix_ctrl_offset, rotation_center_4x4xN)
    if(~exist('fix_ctrl_offset','var')), fix_ctrl_offset = false; end
    
    P = size(p,2);          % number of population size
    N = num_transforms*6*P; 
    M = num_ctrl_pnts;     % number of control points

    p_ctrl_NxM = reshape(permute(reshape(p(1:num_transforms*M*6,:),num_transforms*6,M,P),[1 3 2]),[],M); % rand(N,M)*2-1;
    if(numel(ctrl_scale)==6)
        ctrl_scale_Nx1 = repmat(ctrl_scale(:), num_transforms*P, 1);
    else
        % in case ctrl_scale is specified per transformation
        ctrl_scale_Nx1 = repmat(ctrl_scale(:), P, 1);
    end
    if(fix_ctrl_offset)
        p_offset_Nx1 = zeros(num_transforms*6*P, 1);
        offset_scale_Nx1 = zeros(num_transforms*6*P, 1);
    else
        p_offset_Nx1 = reshape(p((1:num_transforms*6) + num_transforms*M*6,:),[],1); % rand(N,1)*2-1;
        offset_scale_Nx1 = repmat(offset_scale', num_transforms*P, 1);
    end
%     global_NxM = repmat(global_transform_6xM, num_transforms*P, 1);
    
    if(num_parallel==1)
        [trans_param_NxF_array, scaled_param] = spline_interp_param(p_ctrl_NxM, p_offset_Nx1, ctrl_scale_Nx1, offset_scale_Nx1, [], num_frames);
    else
        trans_param_NxF_cell = cell(num_parallel, 1);
        indx_array = array_split_1D(1:N, num_parallel);
        p_ctrl_NxM_par = mat2cell(p_ctrl_NxM, cellfun(@length, indx_array), M);
        p_offset_Nx1_par = mat2cell(p_offset_Nx1, cellfun(@length, indx_array), 1);
        ctrl_scale_Nx1_par = mat2cell(ctrl_scale_Nx1, cellfun(@length, indx_array), 1);
        offset_scale_Nx1_par = mat2cell(offset_scale_Nx1, cellfun(@length, indx_array), 1);
%         global_NxM_par = mat2cell(global_NxM, cellfun(@length, indx_array), M);
        parfor i=1:num_parallel
%        for i=1:num_parallel
            trans_param_NxF_cell{i} = spline_interp_param(p_ctrl_NxM_par{i}, p_offset_Nx1_par{i}, ...
                                            ctrl_scale_Nx1_par{i}, offset_scale_Nx1_par{i}, [], num_frames);
        end
        trans_param_NxF_array = cell2mat(trans_param_NxF_cell);
        scaled_param = [];
    end
    
    if(exist('volume_offset_4x4xN','var'))
        t_6xN = reshape(trans_param_NxF_array,6,[]);
        trans_param_4x4_image_coord = zeros(4, 4, size(t_6xN,2));
        trans_param_4x4_local_coord = zeros(4, 4, size(t_6xN,2));
%         if(isempty(rotation_center_4x4xN))
%             rotation_center_4x4xN = repmat(eye(4), [1 1 num_transforms]);
%         end
        for i=1:num_transforms
            [trans_param_4x4_image_coord(:,:,i:num_transforms:end), trans_param_4x4_local_coord(:,:,i:num_transforms:end)] = ...
                convertTransRotTo4x4_multi_PrePostOffset(t_6xN(:,i:num_transforms:end),volume_offset_4x4xN(:,:,i),inv(volume_offset_4x4xN(:,:,i)),rotation_center_4x4xN(:,:,i));
        end
        trans_param_4x4_image_coord = reshape(trans_param_4x4_image_coord,4,4,num_transforms,P,[]);
        trans_param_4x4_local_coord = reshape(trans_param_4x4_local_coord,4,4,num_transforms,P,[]);
    else
        trans_param_4x4_image_coord = reshape(RegTools.convertTransRotTo4x4_multi(reshape(trans_param_NxF_array,6,[])),4,4,num_transforms,P,[]);
    end
    
    % add global transformation
    for j=1:num_transforms
        global_transform_4x4xM = RegTools.convertTransRotTo4x4_multi(global_transform_6xMxN(:,:,j));
        for i=1:size(global_transform_6xMxN,2)
            trans_param_4x4_image_coord(:,:,j,:,i) = reshape(global_transform_4x4xM(:,:,i)*reshape(trans_param_4x4_image_coord(:,:,j,:,i),4,[]),4,4,1,P);
            trans_param_4x4_local_coord(:,:,j,:,i) = reshape(global_transform_4x4xM(:,:,i)*reshape(trans_param_4x4_local_coord(:,:,j,:,i),4,[]),4,4,1,P);
        end
    end   