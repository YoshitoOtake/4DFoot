function transformation_curvature = compute_transformation_curvature(trans_4x4xNxPxF)

[~, ~, N, P, F] = size(trans_4x4xNxPxF);
if(F<=2), transformation_curvature = 0; return; end

% calculate relative transformation from the previous frame
% (we assume the relative transformation becomes small)
relative_4x4xNxPxF_1 = repmat(eye(4), [1 1 N P F-1]);
for k=1:N
    for l=1:P
        for i=2:F
            relative_4x4xNxPxF_1(:,:,k,l,i-1) = trans_4x4xNxPxF(:,:,k,l,i-1) \ trans_4x4xNxPxF(:,:,k,l,i);
        end
    end
end

% convert relative transformation matrices to 6x1 vectors
trans_param_6xNxPxF_1_array_reshaped = reshape(RegTools.convert4x4ToTransRot_multi(reshape(relative_4x4xNxPxF_1,4,4,[])),6,N, P, F-1);

% difference of the relative transformation, i.e., curvature
euclidean_diff = sqrt(sum( (trans_param_6xNxPxF_1_array_reshaped(:,:,:,2:end)-trans_param_6xNxPxF_1_array_reshaped(:,:,:,1:end-1)).^2,1));
transformation_curvature = reshape(mean(mean(euclidean_diff, 4),2),1,P);
