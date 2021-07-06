function [landmark_tbls, joint_center_tbl, landmark_2D_filename, landmark_3D_filename] = LoadLandmarkFiles_foot(root_dir, x_ray_landmark_dir, patient_ID, view_direction, xray_image_files, x_ray_landmark_filename, bone_label, leading_view, CT_file_name, CT_landmark_filename)
landmark_tbls = cell(2,1);
z_offset = 0; %-(1806-700)/2*0.3;
landmark_2D_filename = cell(2,1);
landmark_3D_filename = cell(2,1);
for i=1:length(view_direction)
    landmark_2D_filename{i} = fullfile(root_dir, 'landmark', x_ray_landmark_dir, patient_ID, view_direction{i}, ['BP_' xray_image_files{i}], x_ray_landmark_filename);
    landmark_3D_filename{i} = fullfile(root_dir, 'landmark', 'CT_left_leg', patient_ID, CT_file_name{1}, view_direction{i}, CT_landmark_filename);
    [pos2D, ~] = LoadDeepLabCutCsvFile(landmark_2D_filename{i});
    pos2D = cleanup_2D_landmark(pos2D); % interpolate NaN
    [pos3D, IDs3D] = LoadSlicerMarkupFile(landmark_3D_filename{i});
    num_frames = size(pos2D,3);
    pos3D = pos3D * [-1 0 0; 0 -1 0; 0 0 1];
    pos3D(:,3) = pos3D(:,3) + z_offset;
%     pos3D = pos3D - ones(size(pos3D,1),1) * (HU_header.Offset.*[-1 -1 1]);  % the origin is at the center of CT bounding box
%     pos3D = pos3D - ones(size(pos3D,1),1) * (HU_header.ElementSpacing.*HU_header.DimSize/2);  % move the origin to center of the CT bounding box
    bone_name = cellfun(@(x) x(1:end-1), IDs3D, 'UniformOutput', false);
    pos2D = cellfun(@(x) squeeze(x)',mat2cell(pos2D, ones(size(pos2D,1),1), size(pos2D,2), size(pos2D,3)),'UniformOutput',false);
    pos2D_cell = cell(size(bone_label,1),1);
    pos3D_cell = cell(size(bone_label,1),1);
    for j=1:size(bone_label,1)
        indx = find(contains(bone_name, bone_label{j}));
        if(isempty(indx))
            pos2D_cell{j} = zeros(num_frames,2,0);  % dummy array
        else
            pos2D_cell{j} = reshape(cell2mat(pos2D(indx)'),[],2,length(indx));
        end
        pos3D_cell{j} = pos3D(indx,:);
    end
    landmark_tbls{i} = table(bone_label(:,1), pos3D_cell, pos2D_cell , 'VariableNames', {'bone_name', 'pos3D', 'pos2D'});
end
% synchronization (one view is delayed by 1 frame)
if(leading_view==0)
    % no modification
elseif(leading_view==1)
    landmark_tbls{1}.pos2D = cellfun(@(x) x(2:end,:,:),landmark_tbls{1}.pos2D,'UniformOutput',false);
    landmark_tbls{2}.pos2D = cellfun(@(x) x(1:end-1,:,:),landmark_tbls{2}.pos2D,'UniformOutput',false);
else
    landmark_tbls{1}.pos2D = cellfun(@(x) x(1:end-1,:,:),landmark_tbls{1}.pos2D,'UniformOutput',false);
    landmark_tbls{2}.pos2D = cellfun(@(x) x(2:end,:,:),landmark_tbls{2}.pos2D,'UniformOutput',false);
end
num_frames = min(size(landmark_tbls{1}.pos2D{1},1), size(landmark_tbls{2}.pos2D{1},1));
landmark_tbls{1}.pos2D = cellfun(@(x) x(1:num_frames,:,:),landmark_tbls{1}.pos2D,'UniformOutput',false);
landmark_tbls{2}.pos2D = cellfun(@(x) x(1:num_frames,:,:),landmark_tbls{2}.pos2D,'UniformOutput',false);

if(1)
    % load joint center landmark
    [pos3D_Joint, IDs3D_Joint] = LoadSlicerMarkupFile( fullfile(root_dir, 'landmark', 'CT_left_leg', patient_ID, CT_file_name{1}, 'joint_centers.fcsv') );
    pos3D_Joint = pos3D_Joint * [-1 0 0; 0 -1 0; 0 0 1];
    pos3D_Joint(:,3) = pos3D_Joint(:,3) + z_offset;
    joint_center_tbl = table(pos3D_Joint, 'VariableNames', {'pos3D'}, 'RowNames', IDs3D_Joint);
else
    local_coordinate = jsondecode( fileread(fullfile(root_dir, 'local_coordinate', patient_ID, 'local_coordinate.json')) );
    joint_center_tbl = table(local_coordinate.center2local_4x4xN(1:3,4,2)'*diag([-1 -1 -1]), 'VariableNames', {'pos3D'}, 'RowNames', {'tibia'});
end