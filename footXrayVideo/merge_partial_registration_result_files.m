addpath(getenv('RegToolsPath'));
addpath('../util');
dataset_ID_array = {'auto', 'manual', 'seg_auto_landmark_manual', 'seg_manual_landmark_auto', 'bone_model_manual', 'bone_model_RSA'};
for dataset_indx = 2 %[1 2 3 4]
    dataset_ID = dataset_ID_array{dataset_indx};
    registration_results_root = fullfile('D:\Collaboration\Nara Medical University\\Foot2D3D_registration_results\', sprintf('20210305_rigidity0_joint0_alpha05_ver1_Foot2D3D_registration_results_%s',dataset_ID)) ;
    output_folder = fullfile(registration_results_root, 'merged');
    if(~exist(output_folder,'dir')), mkdir(output_folder); end
    bone_mode_array = { 'tarsal_bones_proximal'}; %{ 'tarsal_bones_proximal', 'tarsal_bones_distal', 'metatarsal_bones' };
    for patient_indx = [1 2 3 4 5]
    % for patient_indx = 4
        [patient_ID, start_frame, end_frame, all_xray_image_files, calibration_patient_ID] = patient_specific_setup(patient_indx);
        for x_ray_image_indx = 1:length(all_xray_image_files)
            x_ray_image_files = all_xray_image_files{x_ray_image_indx};
            if(patient_indx<0), 
                bone_mode_indx_try = 1; 
        %         patient_ID = 'MODEL1'; 
            else,
                bone_mode_indx_try = 1:length(bone_mode_array); 
        %         patient_ID = sprintf('ND%04d',patient_indx);
            end
            unique_frame_cell = cell(length(bone_mode_indx_try),1);
            num_transforms_array = zeros(length(bone_mode_indx_try),1);
            trans_param_4x4_cell = cell(length(bone_mode_indx_try),1);
            all_max_similarity_cell = cell(length(bone_mode_indx_try),1);
            bone_tbl_cell =cell(length(bone_mode_indx_try),1);

            for k = bone_mode_indx_try
                bone_mode = bone_mode_array{k};
                bone_tbl_cell{k} = generate_bone_table( bone_mode );
                if(patient_indx>0)
                    registration_results_folder = getLatestTimestampedFile( fullfile(registration_results_root, sprintf('results_local_ND%04d_%s_%s_%s*',patient_indx,x_ray_image_files{1},x_ray_image_files{2},bone_mode)) );
                else
                    registration_results_folder = getLatestTimestampedFile( fullfile(registration_results_root, 'results_*local*_MODEL1_*') );
                end
                fprintf('merging result files of %s of ND%04d, in %s\n', bone_mode, patient_indx, registration_results_folder.name);

                all_registration_result_files = dir( fullfile(registration_results_folder.folder, registration_results_folder.name, 'final_results_*.json') );

                all_registration_results_tbl = table([], {}, [], 'VariableNames', {'frame', 'transforms', 'final_similarity_per_frame'});
                for i=1:length(all_registration_result_files)
                    registration_result_file = fullfile( all_registration_result_files(i).folder, all_registration_result_files(i).name );
                    loaded_result = jsondecode(fileread(registration_result_file));
                    loaded_result.transformation_parameters = reshape(loaded_result.transformation_parameters,4,4,loaded_result.num_transforms,[]);
                    if(i==1), volume_offset_4x4xN = reshape(loaded_result.volume_offset_4x4xN,4,4,loaded_result.num_transforms,[]); end
                    for j=1:length(loaded_result.registration_frame)
                        all_registration_results_tbl(end+1,:) = {loaded_result.registration_frame(j), ...
                            RegTools.convert4x4ToTransRot_multi( squeeze(loaded_result.transformation_parameters(:,:,:,j)) ), loaded_result.final_similarity_per_frame(j)};
                    end
                end

                unique_frame_cell{k} = unique(all_registration_results_tbl.frame);
                [num_param, num_transforms_array(k)] = size(all_registration_results_tbl.transforms{1});
                trans_param_4x4_cell{k} = zeros(4,4,num_transforms_array(k),length(unique_frame_cell{k}));
                all_max_similarity_cell{k} = zeros(length(unique_frame_cell{k}),1);
                for i=1:length(unique_frame_cell{k})
                    % pick the frame with maximum similarity when multiple registration results exist
                    all_trans = all_registration_results_tbl.transforms(all_registration_results_tbl.frame==unique_frame_cell{k}(i));
                    all_similarity = all_registration_results_tbl.final_similarity_per_frame(all_registration_results_tbl.frame==unique_frame_cell{k}(i));
                    [all_max_similarity_cell{k}(i), max_indx] = max(all_similarity);
                    trans_param_4x4_cell{k}(:,:,:,i) = RegTools.convertTransRotTo4x4_multi( all_trans{max_indx} );
                end
            end

            unique_frame = unique(cell2mat(unique_frame_cell));
            bone_names = cellfun(@(x) x.Row', bone_tbl_cell, 'UniformOutput', false);
            bone_names = [bone_names{:}];
            bone_indx = mat2cell(1:length(bone_names), 1, cellfun(@(x) length(x.Row), bone_tbl_cell)');
            bone_label_indx = cell(length(bone_names),1);
            trans_param_4x4 = zeros(4, 4, length(bone_names), length(unique_frame));
            all_max_similarity = zeros(length(unique_frame), length(bone_names));
            for j=1:length(trans_param_4x4_cell)
                trans_param_4x4(:,:,bone_indx{j},:) = trans_param_4x4_cell{j};
                all_max_similarity(:,bone_indx{j}) = repmat(all_max_similarity_cell{j},[1 length(bone_indx{j})]);
                bone_label_indx(bone_indx{j}) = bone_tbl_cell{j}.label_indx;
            end

            if(0)
                %%
                figure('Position',[100 100 [1600 900]*8/10], 'PaperPositionMode', 'auto', 'Color', 'w', 'InvertHardcopy', 'off');
                similarity = cell2mat(all_max_similarity_cell');
                hold on;
                plot(similarity);
                sample_X = floor(linspace(1,size(similarity,1),5));
                plot(interp1(sample_X,similarity(sample_X,:),1:size(similarity,1),'spline'));

                %%
                relative_4x4 = repmat(eye(4), [1 1 size(trans_param_4x4,3) size(trans_param_4x4,4)]);
                for i=1:length(bone_names)
                    for j=2:length(unique_frame)
                        relative_4x4(:,:,i,j) = inv(trans_param_4x4(:,:,i,1))*trans_param_4x4(:,:,i,j);
                    end
                end

               relative_6xN = reshape( RegTools.convert4x4ToTransRot_multi( reshape(relative_4x4, 4, 4, []) ), 6, length(bone_names), length(unique_frame));
                figure('Position',[100 100 [1600 900]*8/10], 'PaperPositionMode', 'auto', 'Color', 'w', 'InvertHardcopy', 'off');
                for i=1:length(bone_names)
                    subplot(3,5,i);
                    plot(squeeze(relative_6xN(:,i,:))');
                    title(sprintf('%s', bone_names{i}), 'interpreter', 'none');
                end
            end
    %%

            selected_bones = 1:length(bone_names);
            S = struct('registration_frame', unique_frame, 'num_transforms', length(bone_names(selected_bones)), 'num_frames', length(unique_frame), ...
               'transformation_parameters', trans_param_4x4(:,:,selected_bones,:), 'final_similarity_per_frame', all_max_similarity(:,selected_bones), ...
               'bone_names', {bone_names(selected_bones)'}, 'bone_label_indx', {bone_label_indx(selected_bones)}, 'volume_offset_4x4xN', volume_offset_4x4xN);

            temp = strcat(bone_mode_array(bone_mode_indx_try), '_');
            fid = fopen( fullfile(output_folder, sprintf('%s_%s_%s_%s', patient_ID, x_ray_image_files{1}, x_ray_image_files{2}, [[temp{:}] 'final_results.json'])), 'w' );
            fprintf(fid, '%s', prettyjson(jsonencode(S)));
            fclose(fid);
        end
    end
end