function [matrix_4x4xN_pre_post_offset, matrix_4x4xN_pre_offset] = convertTransRotTo4x4_multi_PrePostOffset(p_6xN, offset_4x4_pre, offset_4x4_post, rotation_center_local)

N = size(p_6xN,2);
zeros_1xN = zeros(1,N);
c = cos(pi/180.0*p_6xN(4:6,:)); s = sin(pi/180.0*p_6xN(4:6,:));
% multiply offset from left-hand side
matrix_16xN = [ c(2,:).*c(3,:);                         c(2,:).*s(3,:); -s(2,:); zeros_1xN;
      -c(1,:).*s(3,:)+s(1,:).*s(2,:).*c(3,:); c(1,:).*c(3,:)+s(1,:).*s(2,:).*s(3,:); s(1,:).*c(2,:); zeros_1xN;
      s(1,:).*s(3,:)+c(1,:).*s(2,:).*c(3,:); -s(1,:).*c(3,:)+c(1,:).*s(2,:).*s(3,:); c(1,:).*c(2,:); zeros_1xN;
      p_6xN(1,:);                            p_6xN(2,:);                            p_6xN(3,:);      ones(1,N)];
temp_4x4xN = reshape( offset_4x4_pre*rotation_center_local*reshape( matrix_16xN, 4,4*N ), 4,4,N);
matrix_4Nx4_pre_offset = reshape(permute(temp_4x4xN,[1 3 2]),4*N,4) / rotation_center_local;
matrix_4x4xN_pre_offset = permute( reshape( matrix_4Nx4_pre_offset, 4, N, 4), [1 3 2]);
matrix_4x4xN_pre_post_offset = permute( reshape( matrix_4Nx4_pre_offset*offset_4x4_post, 4, N, 4), [1 3 2]);
