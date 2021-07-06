function [trans_param_NxF, scaled_param] = spline_interp_param(p_ctrl_NxM, p_offset_Nx1, ctrl_scale_Nx1, offset_scale_Nx1, global_NxM, num_frames)

N = size(p_ctrl_NxM, 1);    % number of parameters to interpolate
M = size(p_ctrl_NxM, 2);    % number of control points

% scaled_param = global_NxM + repmat(p_offset_Nx1.*offset_scale_Nx1, 1, M) + diag(ctrl_scale_Nx1)*p_ctrl_NxM;
scaled_param = repmat(p_offset_Nx1.*offset_scale_Nx1, 1, M) + diag(ctrl_scale_Nx1)*p_ctrl_NxM;
if(M==num_frames)
    trans_param_NxF = scaled_param;
else
    x = linspace(0,1,M);
    xq = linspace(0,1,num_frames);
    F = griddedInterpolant(x, 1:M, 'spline');
    trans_param_NxF = zeros(N, num_frames);
    for i=1:N
    %     trans_param_NxF(i,:) = interp1(x, scaled_param(i,:), xq, 'spline');
        F.Values = scaled_param(i,:);
        trans_param_NxF(i,:) = F(xq);
    end
end
