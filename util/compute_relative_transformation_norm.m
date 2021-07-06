function relative_transformation_norm = compute_relative_transformation_norm(trans_4x4xNxPxF)

[~, ~, N, P, F] = size(trans_4x4xNxPxF);

% calculate relative transformation with respect to the first object
% (we assume the relative transformation becomes small)
relative_4x4xN_1xPxF = repmat(eye(4), [1 1 N-1 P F]);
for k=1:P
    for l=1:F
        for i=2:N
            relative_4x4xN_1xPxF(:,:,i-1,k,l) = trans_4x4xNxPxF(:,:,1,k,l) \ trans_4x4xNxPxF(:,:,i,k,l);
        end
    end
end

% convert relative transformation matrices to 6x1 vectors
relative_6xN_1xFxP_array_reshaped = reshape(RegTools.convert4x4ToTransRot_multi(reshape(relative_4x4xN_1xPxF,4,4,[])),6,N-1, P, F);

% variance with respect to the identity matrix, (0,0,0,0,0,0)
relative_transformation_norm = reshape(mean(mean(sum(relative_6xN_1xFxP_array_reshaped.^2,1),2),4),1,P);
