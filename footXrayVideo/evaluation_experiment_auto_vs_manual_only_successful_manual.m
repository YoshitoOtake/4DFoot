addpath(getenv('RegToolsPath'));
results_root_dir = 'D:\Collaboration\Nara Medical University\Foot2D3D_registration_results\';
% date_str = '20210220_ver1';
% date_str = '20210221_rigidity5e-4';
date_str = '20210302_rigidity1e-3_alpha1';
data_folder_prefix = [date_str '_Foot2D3D_registration_results'];
dataset_ID_array = {'auto', 'manual', 'seg_auto_landmark_manual', 'seg_manual_landmark_auto', 'bone_model_manual', 'bone_model_RSA'};
dataset_indx_array = 2; %[1 3 4];
all_stats = zeros(length(dataset_indx_array), 12);
for dataset_indx = 1:length(dataset_indx_array)
    dataset_ID = dataset_ID_array{dataset_indx_array(dataset_indx)};
    auto_dir = fullfile( results_root_dir, [data_folder_prefix, '_', dataset_ID], 'merged');
%     auto_dir = fullfile( results_root_dir, ['20210221_rigidity5e-4_Foot2D3D_registration_results', '_', dataset_ID], 'merged');
%     manual_dir = fullfile( results_root_dir, [data_folder_prefix, '_manual'], 'merged');
    manual_dir = fullfile( results_root_dir, '20210220_ver1_Foot2D3D_registration_results_manual', 'merged');
%     manual_summary_file_prefix = fullfile( results_root_dir, ['summary_' date_str], [data_folder_prefix, '_manual_success_summary'] );
   manual_summary_file_prefix = fullfile( results_root_dir, 'summary_20210219', '20210219_ver1_Foot2D3D_registration_results_manual_success_summary' );
   output_folder = fullfile( results_root_dir, ['summary_' date_str], 'all_patients_in_one_plot');
    if(~exist(output_folder,'dir')), mkdir(output_folder); end

    % bone_list = jsondecode(fileread('ankle_new.json'));
    % bone_names = cellfun(@(x) x{2},bone_list.color,'UniformOutput',false);
    % bone_label = bone_names([1 3:14]); %{'tibia', 'talus', 'calcaneal', 'navicular', 'Medial_cuneiform', 'Intermediate_cuneiform', 'Lateral_cuneiform', 'cuboid', '1st_metatarsal', '}';
    bone_label = num2cell(1:13)';
    num_bones = length(bone_label);

    patient_indx_array = [1 2 3 4 5];
    all_error_6xN_cell = cell(length(patient_indx_array),1);
    fig1 = figure('Position',[300 300 [1200 800]*8/10], 'PaperPositionMode', 'auto', 'Color', 'w');
    all_summary_linear_success_only = cell(3, length(patient_indx_array));
    all_success_flags = cell(1, length(patient_indx_array));
    for patient_indx = patient_indx_array
        [patient_ID, start_frame, end_frame, all_xray_image_files, calibration_patient_ID] = patient_specific_setup(patient_indx);
        summary_tbl = readtable( [manual_summary_file_prefix, '_', patient_ID, '.csv'] );
        all_success_flags{patient_indx} = summary_tbl{:,2:end};
        registration_frames = summary_tbl{:,1};
        num_successful_frames = sum(all_success_flags{patient_indx});

        local_coordinate_file = fullfile('\\scallop\User\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data\local_coordinate\', patient_ID, 'local_coordinate.json');
        local_coordinate = jsondecode( fileread(local_coordinate_file) );
    %     for j=1:size(local_coordinate.center2local_4x4xN,3), local_coordinate.center2local_4x4xN(1:3,4,j) = local_coordinate.center2local_4x4xN(1:3,4,j) * -1; end    % wrt local coordinate
        for j=1:size(local_coordinate.center2local_4x4xN,3), local_coordinate.center2local_4x4xN(:,:,j) = inv(local_coordinate.center2local_4x4xN(:,:,j)) * diag([-1 -1 1 1]); end    % wrt local coordinate

        auto_json = dir( fullfile(auto_dir, sprintf('ND%04d*.json',patient_indx)) );
        manual_json = dir( fullfile(manual_dir,sprintf('ND%04d*.json',patient_indx)) );

        auto_data = jsondecode( fileread(fullfile(auto_json(1).folder, auto_json(1).name)) );
        manual_data = jsondecode( fileread(fullfile(manual_json(1).folder, manual_json(1).name)) );

        all_error_6xN = zeros([length(registration_frames), size(all_success_flags{patient_indx},2), 6]);
        all_error_translation = zeros([length(registration_frames), size(all_success_flags{patient_indx},2)]);
        all_error_rotation = zeros([length(registration_frames), size(all_success_flags{patient_indx},2)]);
        for i=1:length(registration_frames)
            for j=1:size(all_success_flags{patient_indx},2)
                if(j<=size(local_coordinate.center2local_4x4xN,3)), local_coord = local_coordinate.center2local_4x4xN(:,:,j); else, local_coord = eye(4); end
                error_trans = inv(manual_data.transformation_parameters(:,:,j,i)*local_coord)*auto_data.transformation_parameters(:,:,j,i)*local_coord;
                all_error_6xN(i,j,:) = RegTools.convert4x4ToTransRot_multi( error_trans );
                all_error_translation(i,j) = norm(error_trans(1:3,4));
                all_error_rotation(i,j) = rad2deg(norm(RotMat2RodAng(error_trans(1:3,1:3))));
            end
        end

        all_summary_linear = {reshape(repmat(bone_label', length(registration_frames), 1), [], 1), all_error_translation(:), all_error_rotation(:)};
        all_summary_linear_success_only(:,patient_indx) = all_summary_linear;
        for i=1:3
            all_summary_linear_success_only{i,patient_indx}(~all_success_flags{patient_indx}(:)) = [];
        end

        %%
        [hs,vs,tb,bb,lb,rb] = deal(0.05,0.05,0.05,0.05,0.05,0.05);

        trans_plot_range = [0 40];
        rotation_plot_range = [0 30];
        msubplot(length(patient_indx_array),3,(patient_indx-1)*3+1, hs,vs,tb,bb,lb,rb);
        bar(mean(all_success_flags{patient_indx})*100);
        set(gca,'ylim',[0 100]);

        msubplot(length(patient_indx_array),3,(patient_indx-1)*3+2, hs,vs,tb,bb,lb,rb);
        boxplot(all_summary_linear_success_only{2,patient_indx}, all_summary_linear_success_only{1,patient_indx});
        set(gca,'ylim',trans_plot_range);

        msubplot(length(patient_indx_array),3,(patient_indx-1)*3+3, hs,vs,tb,bb,lb,rb);
        boxplot(all_summary_linear_success_only{3,patient_indx}, all_summary_linear_success_only{1,patient_indx});
        set(gca,'ylim',rotation_plot_range);
    end
    annotation('textbox', [0 0.98 1/3 0.02],'String','Success rate (%)','HorizontalAlignment','Center','FontSize',18,'LineStyle','none');
    annotation('textbox', [1/3 0.98 1/3 0.02],'String','Translation error (mm)','HorizontalAlignment','Center','FontSize',18,'LineStyle','none');
    annotation('textbox', [2/3 0.98 1/3 0.02],'String','Rotation error (deg)','HorizontalAlignment','Center','FontSize',18,'LineStyle','none');

    print(gcf, '-dpng', '-r200', fullfile(output_folder, sprintf('ND0001_ND0005_full_auto_registration_error_summary_%s.png', dataset_ID)));

    %%
    if(1)
        figure('Position',[300 300 [800 700]*8/10], 'PaperPositionMode', 'auto', 'Color', 'w');
        all_success_flags_linear = cell2mat(cellfun(@(x) x', all_success_flags, 'UniformOutput', false))';
        bar(mean(all_success_flags_linear)*100);
        set(gca,'ylim',[0 100], 'FontSize', 20);
        ylabel('Success rate (%)');
        print(gcf, '-dpng', '-r200', fullfile(output_folder, sprintf('ND0001_ND0005_full_auto_registration_error_summary_success_rate_in_one_plot_%s.png', dataset_ID)));
    end
    %%
    if(1)
        figure('Position',[300 300 [1400 700]*8/10], 'PaperPositionMode', 'auto', 'Color', 'w');
        [hs,vs,tb,bb,lb,rb] = deal(0.07,0.05,0.08,0.05,0.05,0.05);
        all_labels = cell(length(patient_indx_array),1);
        for patient_indx = 1:length(patient_indx_array)
            all_labels = [all_labels; all_summary_linear_success_only{1,patient_indx}];
        end
        all_labels = cell2mat( all_labels );
        all_translations = cell2mat(all_summary_linear_success_only(2,:)');
        all_rotations = cell2mat(all_summary_linear_success_only(3,:)');
        msubplot(1,2,1, hs,vs,tb,bb,lb,rb);
        boxplot( all_translations, all_labels );
    %     title('Translation (mm), All patients', 'FontSize', 18);
        set(gca,'ylim',[0 8], 'FontSize', 20);
        ylabel('Translation (mm)');
        msubplot(1,2,2, hs,vs,tb,bb,lb,rb);
        boxplot( all_rotations, all_labels );
    %     title('Rotation (deg), All patients', 'FontSize', 18);
        set(gca,'ylim',[0 4], 'FontSize', 20);
        ylabel('Rotation (deg)');
        print(gcf, '-dpng', '-r200', fullfile(output_folder, sprintf('All_patients_full_auto_registration_error_summary_%s.png', dataset_ID)));

%         all_stats = zeros(4, num_bones);
%         for i=1:num_bones
%             all_stats(1,i) = mean( all_translations(all_labels==i) );
%             all_stats(2,i) = std( all_translations(all_labels==i) );
%             all_stats(3,i) = mean( all_rotations(all_labels==i) );
%             all_stats(4,i) = std( all_rotations(all_labels==i) );
%         end
    end
    all_stats(dataset_indx,:) = [mean( all_translations(ismember(all_labels,2:4) ) ), std( all_translations(ismember(all_labels,2:4) ) ), mean( all_rotations(ismember(all_labels,2:4) ) ), std( all_rotations(ismember(all_labels,2:4) ) ), ...
        mean( all_translations(ismember(all_labels,5:8) ) ), std( all_translations(ismember(all_labels,5:8) ) ), mean( all_rotations(ismember(all_labels,5:8) ) ), std( all_rotations(ismember(all_labels,5:8) ) ), ...
        mean( all_translations(ismember(all_labels,9:13) ) ), std( all_translations(ismember(all_labels,9:13)) ), mean( all_rotations(ismember(all_labels,9:13) ) ), std( all_rotations(ismember(all_labels,9:13)) ) ...
    ];
end