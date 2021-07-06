% root = '\\scallop\user\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data\local_coordinate';
% all_local_coordinate = zeros(4,4,4,6);
% for patient_indx=1:5
%         patient_ID = sprintf('ND%04d',patient_indx);
%         local_coordinate = jsondecode( fileread(fullfile(root, patient_ID, 'local_coordinate.json')) );
%         all_local_coordinate(:,:,:,patient_indx) = local_coordinate.center2local_4x4xN;
% end
% patient_ID = 'MODEL1'; 
% root_model = '\\scallop\user\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data_model\local_coordinate';
% local_coordinate = jsondecode( fileread(fullfile(root_model, patient_ID, 'local_coordinate.json')) );
% all_local_coordinate(:,:,:,6) = local_coordinate.center2local_4x4xN;
% 
% return;

Xray_suffix = '';
dataset_ID_array = {'auto', 'manual', 'semi-auto', 'bone_model_manual', 'bone_model_RSA'};
for patient_indx = [1 2 3 4 5]
    for dataset_indx = [2 1]
        dataset_ID = dataset_ID_array{dataset_indx};
        fprintf('start footXrayVideo_NaraMedUniv, dataset_ID: %s\n', dataset_ID);
        setup_dataset;
        patient_ID = sprintf('ND%04d',patient_indx);
        local_coordinate = jsondecode( fileread(fullfile(data_root_dir_landmark, 'local_coordinate', patient_ID, 'local_coordinate.json')) );
        talus_coordinate = inv(local_coordinate.center2local_4x4xN(:,:,2));
        out_filename = fullfile(data_root_dir_landmark, 'landmark', 'CT_left_leg', patient_ID, CT_file_name{1}, 'joint_centers.fcsv');
%         SaveSlicerMarkupFile(out_filename, local_coordinate.center2local_4x4xN(1:3,4,2)'*diag([1 1 -1]), {'tibia'});
        SaveSlicerMarkupFile(out_filename, talus_coordinate(1:3,4)'*diag([-1 -1 1]), {'tibia'});
    end
end

%%
patient_ID = 'MODEL1'; 
dataset_ID = dataset_ID_array{4};
setup_dataset;

local_coordinate = jsondecode( fileread(fullfile(data_root_dir_landmark, 'local_coordinate', patient_ID, 'local_coordinate.json')) );
talus_coordinate = inv(local_coordinate.center2local_4x4xN(:,:,2));
out_filename = fullfile(data_root_dir_landmark, 'landmark', 'CT_left_leg', patient_ID, CT_file_name{1}, 'joint_centers.fcsv');
% SaveSlicerMarkupFile(out_filename, local_coordinate.center2local_4x4xN(1:3,4,2)'*diag([1 1 -1]), {'tibia'});
SaveSlicerMarkupFile(out_filename,talus_coordinate(1:3,4)'*diag([-1 -1 1]), {'tibia'});
