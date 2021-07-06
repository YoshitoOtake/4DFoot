function global_transformation_6xMxN = getGlobalTransformation_from_previous_registration( registration_results_root_dir, patient_ID, xray_image_files, bone_mode, ref_bone_indx, num_total_frames, num_bones)

local_registration_result_dir = getLatestTimestampedFile( fullfile(registration_results_root_dir, sprintf('results_local_%s_%s_%s_%s_*',patient_ID,xray_image_files{1},xray_image_files{2},bone_mode)) );
all_registration_result_files = dir( fullfile(local_registration_result_dir.folder, local_registration_result_dir.name, 'final_results_*.json') );

all_transformation = repmat(eye(4), [1, 1, num_total_frames]);
for i=1:length(all_registration_result_files)
    registration_result_file = fullfile( all_registration_result_files(i).folder, all_registration_result_files(i).name );
    loaded_result = jsondecode(fileread(registration_result_file));
    loaded_result.transformation_parameters = reshape(loaded_result.transformation_parameters,4,4,loaded_result.num_transforms,[]);
    all_transformation(:,:,loaded_result.registration_frame) = loaded_result.transformation_parameters(:,:,ref_bone_indx);
end

global_transformation_6xMxN = repmat( reshape( RegTools.convert4x4ToTransRot_multi(all_transformation), 6, [] ), [1 1 num_bones]);
