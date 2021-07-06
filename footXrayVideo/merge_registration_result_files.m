addpath(getenv('RegToolsPath'));
addpath('../util');
% registration_results_folder = 'D:\Collaboration\Nara Medical University\20200910_Foot2D3D\Matlab_script/results_ND0001_LL0028_LL0009_20201027_221600_success';
% registration_results_root = 'D:\Collaboration\Nara Medical University\20210116_Foot2D3D_registration_results\20210125_MODEL1_success_manual_landmarks';
% registration_results_root = 'D:\Collaboration\Nara Medical University\20210116_Foot2D3D_registration_results\20210127_MODEL1_success_RSA';
registration_results_root = 'D:\Collaboration\Nara Medical University\20210214_ver1_Foot2D3D_registration_results_auto';
output_folder = fullfile(registration_results_root, 'merged');
if(~exist(output_folder,'dir')), mkdir(output_folder); end
for patient_indx = 2 %-1 %4
    if(patient_indx>0)
        registration_results_folder = getLatestTimestampedFile( fullfile(registration_results_root, sprintf('results_*_ND%04d_*',patient_indx)) );
    else
        registration_results_folder = getLatestTimestampedFile( fullfile(registration_results_root, 'results_*local*_MODEL1_*') );
    end
    fprintf('merging result files of ND%04d, in %s\n', patient_indx, registration_results_folder.name);

    all_registration_result_files = dir( fullfile(registration_results_folder.folder, registration_results_folder.name, 'final_results_*.json') );

    all_registration_results_tbl = table([], {}, [], 'VariableNames', {'frame', 'transforms', 'final_similarity_per_frame'});
    for i=1:length(all_registration_result_files)
        registration_result_file = fullfile( all_registration_result_files(i).folder, all_registration_result_files(i).name );
        loaded_result = jsondecode(fileread(registration_result_file));
        loaded_result.transformation_parameters = reshape(loaded_result.transformation_parameters,4,4,loaded_result.num_transforms,[]);
        for j=1:length(loaded_result.registration_frame)
            all_registration_results_tbl(end+1,:) = {loaded_result.registration_frame(j), ...
                RegTools.convert4x4ToTransRot_multi( squeeze(loaded_result.transformation_parameters(:,:,:,j)) ), loaded_result.final_similarity_per_frame(j)};
        end
    end

    unique_frame = unique(all_registration_results_tbl.frame);
    [num_param, num_transforms] = size(all_registration_results_tbl.transforms{1});
    trans_param_4x4 = zeros(4,4,num_transforms,length(unique_frame));
    all_max_similarity = zeros(length(unique_frame),1);
    for i=1:length(unique_frame)
        % pick the frame with maximum similarity when multiple registration results exist
        all_trans = all_registration_results_tbl.transforms(all_registration_results_tbl.frame==unique_frame(i));
        all_similarity = all_registration_results_tbl.final_similarity_per_frame(all_registration_results_tbl.frame==unique_frame(i));
        [all_max_similarity(i), max_indx] = max(all_similarity);
        trans_param_4x4(:,:,:,i) = RegTools.convertTransRotTo4x4_multi( all_trans{max_indx} );
%         for j=1:length(all_trans)
%             rot = all_trans{j}(4:6,:);
%             rot(rot<0) = rot(rot<0)+360;
%             rot(rot>360) = rot(rot>360)-360;
%             all_trans{j}(4:6,:) = rot;
%         end
%         trans_param_4x4(:,:,:,unique_frame(i)) = RegTools.convertTransRotTo4x4_multi( all_trans{ceil(length(all_trans)/2)} );
    %     trans_param_4x4(:,:,:,unique_frame(i)) = RegTools.convertTransRotTo4x4_multi( reshape(mean(cell2mat(cellfun(@(x) reshape(x,1,[]),all_trans,'UniformOutput',false)),1),num_param, num_transforms) );
    end

    S = struct('registration_frame', unique_frame, 'num_transforms', num_transforms, 'num_frames', length(unique_frame), ...
       'transformation_parameters', trans_param_4x4, 'final_similarity_per_frame', all_max_similarity);

    fid = fopen( fullfile(output_folder, [registration_results_folder.name '_final_results.json']), 'w' );
    fprintf(fid, '%s', prettyjson(jsonencode(S)));
    fclose(fid);
end