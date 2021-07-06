addpath(getenv('RegToolsPath'));
results_root_dir = 'D:\Collaboration\Nara Medical University\Foot2D3D_registration_results\';
date_str = '20210305_rigidity0_joint0_alpha05_ver1';
output_folder = fullfile(results_root_dir, ['summary_' date_str]);
if(~exist(output_folder,'dir')), mkdir(output_folder); end
data_folder_prefix = [date_str '_Foot2D3D_registration_results'];
dataset_ID_array = {'auto', 'manual', 'seg_auto_landmark_manual', 'seg_manual_landmark_auto', 'bone_model_manual', 'bone_model_RSA'};
dataset_indx_array = 2; %[1 3 4];
all_stats = zeros(length(dataset_indx_array), 12);
for dataset_indx = 1:length(dataset_indx_array)
    dataset_ID = dataset_ID_array{dataset_indx_array(dataset_indx)};
    auto_dir = fullfile( results_root_dir, [data_folder_prefix, '_', dataset_ID], 'merged');
%     auto_dir = fullfile( results_root_dir, ['20210221_rigidity5e-4_Foot2D3D_registration_results', '_', dataset_ID], 'merged');
%     manual_dir = fullfile( results_root_dir, [data_folder_prefix, '_manual'], 'merged');
%     manual_dir = fullfile( results_root_dir, '20210220_ver1_Foot2D3D_registration_results_manual', 'merged');
%     manual_summary_file_prefix = fullfile( results_root_dir, ['summary_' date_str], [data_folder_prefix, '_manual_success_summary'] );
%    manual_summary_file_prefix = fullfile( results_root_dir, 'summary_20210219', '20210219_ver1_Foot2D3D_registration_results_manual_success_summary' );

    % bone_list = jsondecode(fileread('ankle_new.json'));
    % bone_names = cellfun(@(x) x{2},bone_list.color,'UniformOutput',false);
    % bone_label = bone_names([1 3:14]); %{'tibia', 'talus', 'calcaneal', 'navicular', 'Medial_cuneiform', 'Intermediate_cuneiform', 'Lateral_cuneiform', 'cuboid', '1st_metatarsal', '}';
%     bone_label = num2cell(1:4)';
%     num_bones = length(bone_label);

    patient_indx_array = [1 2 3 4 5];
    for patient_indx = patient_indx_array
        [patient_ID, start_frame, end_frame, all_xray_image_files, calibration_patient_ID] = patient_specific_setup(patient_indx);
        for x_ray_image_indx = 1:length(all_xray_image_files)
            x_ray_image_files = all_xray_image_files{x_ray_image_indx};

            local_coordinate_file = fullfile('\\scallop\User\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data\local_coordinate\', patient_ID, 'local_coordinate.json');
            local_coordinate = jsondecode( fileread(local_coordinate_file) );
%             for j=1:size(local_coordinate.center2local_4x4xN,3), local_coordinate.center2local_4x4xN(:,:,j) = inv(local_coordinate.center2local_4x4xN(:,:,j)) * diag([-1 -1 1 1]); end    % wrt local coordinate

            auto_json = dir( fullfile(auto_dir, sprintf('ND%04d_%s_%s_*.json',patient_indx,x_ray_image_files{1},x_ray_image_files{2})) );

            auto_data = jsondecode( fileread(fullfile(auto_json(1).folder, auto_json(1).name)) );
            local_coordinates_all_frames = zeros(4,4,size(local_coordinate.center2local_4x4xN,3),length(auto_data.registration_frame));
            for i=1:length(auto_data.registration_frame)
                for j=1:size(local_coordinate.center2local_4x4xN,3)
                    local_coordinates_all_frames(:,:,j,i) = auto_data.transformation_parameters(:,:,j,i)*inv(local_coordinate.center2local_4x4xN(:,:,j));
                end
            end
            joint_list = [1 2; 2 3; 2 4];
            joint_angles = zeros(6, size(joint_list,1), length(auto_data.registration_frame));
            for i=1:length(auto_data.registration_frame)
                for j=1:size(joint_list,1)
                    k = joint_list(j,:);
                    joint_angles(:,j,i) = RegTools.convert4x4ToTransRot_multi( inv(local_coordinates_all_frames(:,:,k(1),i)) * local_coordinates_all_frames(:,:,k(2),i) );
                end
            end
            
            fid = fopen( fullfile(output_folder, sprintf('ND%04d_%s_%s.csv',patient_indx,x_ray_image_files{1},x_ray_image_files{2})), 'w' );
            fprintf(fid, 'frame no,');
            for j=1:size(joint_list,1)
                fprintf(fid, '(joint%d) xtrans,ytrans,ztrans,xrot,yrot,zrot', j);
                if(j==size(joint_list,1))
                    fprintf(fid, '\n');
                else
                    fprintf(fid, ',');
                end
            end

            for i=1:length(auto_data.registration_frame)
                fprintf(fid, '%d,', auto_data.registration_frame(i));
                for j=1:size(joint_list,1)
                    fprintf(fid, '%f, %f, %f, %f, %f, %f', joint_angles(:,j,i));
                    if(j==size(joint_list,1))
                        fprintf(fid, '\n');
                    else
                        fprintf(fid, ',');
                    end
                end
            end
            fclose(fid);
        end
    end
end