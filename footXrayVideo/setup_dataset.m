view_direction = {'lateral', 'oblique'};
if(~exist('enable_SuperSloMo_interpolation','var')), enable_SuperSloMo_interpolation = true; end

if(enable_SuperSloMo_interpolation)
    leading_view = 2; % the "leading_view" starts earlier than the other view by 1 frame
    Xray_suffix = '_interpolated';
else
    leading_view = 0; % no "leading_view", i.e., use the same frame index for both views
    Xray_suffix = '';
end

switch dataset_ID
    case 'auto'
        % 1) auto landmarks
        data_root_dir_image = '\\scallop\user\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data_predicted'; 
        data_root_dir_landmark = '\\scallop\user\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data_predicted'; 
        x_ray_landmark_dir = sprintf('Xray_DLC_csv%s',Xray_suffix); 
        if(isempty(Xray_suffix)), x_ray_landmark_filename = 'dlc_original.csv'; else, x_ray_landmark_filename = 'dlc_interp.csv'; end
        CT_file_name = {'1500', 'image.mhd'}; CT_landmark_filename = 'landmark.fcsv'; x_ray_image_file_suffix = '';
        fixed_clim = [0 1.5];  grad_clim = [-0.005 0.005]; moving_clim = [0 20];
    case 'manual'
        % 2) manual landmarks
        data_root_dir_image = '\\scallop\user\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data'; 
        data_root_dir_landmark = '\\scallop\user\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data'; 
        x_ray_landmark_dir = sprintf('Xray_csv%s',Xray_suffix); x_ray_landmark_filename = 'CollectedData_miyamoto.csv';  
        CT_file_name = {'1500', 'image.mhd'}; CT_landmark_filename = 'landmark.fcsv'; x_ray_image_file_suffix = '';
        fixed_clim = [0 1.5];  grad_clim = [-0.005 0.005]; moving_clim = [0 20];
    case 'seg_auto_landmark_manual'
        % 3) auto segmentation, manual landmark
        data_root_dir_image = '\\scallop\user\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data_predicted'; 
        data_root_dir_landmark = '\\scallop\user\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data'; 
        x_ray_landmark_dir = sprintf('Xray_csv%s',Xray_suffix); x_ray_landmark_filename = 'CollectedData_miyamoto.csv'; 
        CT_file_name = {'1500', 'image.mhd'}; CT_landmark_filename = 'landmark.fcsv'; x_ray_image_file_suffix = '';
        fixed_clim = [0 1.5];  grad_clim = [-0.005 0.005]; moving_clim = [0 20];
    case 'seg_manual_landmark_auto'
        % 4) manual segmentation, auto landmark
        data_root_dir_image = '\\scallop\user\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data'; 
        data_root_dir_landmark = '\\scallop\user\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data_predicted'; 
        x_ray_landmark_dir = sprintf('Xray_DLC_csv%s',Xray_suffix); 
        if(isempty(Xray_suffix)), x_ray_landmark_filename = 'dlc_original.csv'; else, x_ray_landmark_filename = 'dlc_interp.csv'; end
        CT_file_name = {'1500', 'image.mhd'}; CT_landmark_filename = 'landmark.fcsv'; x_ray_image_file_suffix = '';
        fixed_clim = [0 1.5];  grad_clim = [-0.005 0.005]; moving_clim = [0 20];
    case 'bone_model_manual'
        % 5) bone model manual landmarks
        data_root_dir_image = '\\scallop\user\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data_model'; 
        data_root_dir_landmark = '\\scallop\user\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data_model'; 
        x_ray_landmark_dir = sprintf('Xray_csv%s',Xray_suffix); x_ray_landmark_filename = 'CollectedData_miyamoto.csv'; x_ray_image_file_suffix = '_inpainted';
        CT_file_name = {'801', 'image_inpainted.mhd'}; CT_landmark_filename = 'landmark.fcsv';
        fixed_clim = [0.6 1.2]; grad_clim = [-0.001 0.001]; moving_clim = [0 7];
    case 'bone_model_RSA'
        % 6) bone model metallic beads
        data_root_dir_image = '\\scallop\user\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data_model'; 
        data_root_dir_landmark = '\\scallop\user\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data_model'; 
        x_ray_landmark_dir = 'metal_csv_interpolated_ver2'; x_ray_landmark_filename = 'CollectedData_miyamoto.csv'; x_ray_image_file_suffix = '';
        CT_file_name = {'801', 'image.mhd'}; CT_landmark_filename = 'landmark_beads.fcsv';
        fixed_clim = [0.6 1.2]; grad_clim = [-0.001 0.001]; moving_clim = [0 7];
end

