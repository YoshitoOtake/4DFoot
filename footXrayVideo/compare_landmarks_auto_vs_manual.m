
dataset_ID_array = {'auto', 'manual', 'bone_model_manual', 'bone_model_RSA'};
all_landmark_tbls = cell(2,1);
for patient_indx = 1:5
    patient_ID = sprintf('ND%04d',patient_indx);
    output_folder = fullfile('\\scallop\user\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data_predicted\landmark_error\Xray_DLC_csv_interpolated', patient_ID);
    if(~exist(output_folder,'dir')), mkdir(output_folder); end
    
    switch patient_indx
        case 1, all_xray_image_files = { {'LL0027', 'LL0008'} }; calibration_patient_ID = patient_ID;
        case 2, all_xray_image_files = { {'LL0048', 'LL0020'} }; calibration_patient_ID = 'ND0002-5';
        case 3, all_xray_image_files = { {'LL0022', 'LL0006'} }; calibration_patient_ID = 'ND0002-5';
        case 4, all_xray_image_files = { {'LL0023', 'LL0007'} }; calibration_patient_ID = 'ND0002-5';
        case 5, all_xray_image_files = { {'LL0016', 'LL0004'} }; calibration_patient_ID = 'ND0002-5';
    end
    for dataset_indx = [1 2]
        dataset_ID = dataset_ID_array{dataset_indx};
        fprintf('loading landmarks, dataset_ID: %s\n', dataset_ID);
        setup_dataset;
        % load landmark files
        view_direction = {'lateral', 'oblique'};
        leading_view = 2; % the "leading_view" starts earlier than the other view by 1 frame
        bone_tbl = generate_bone_table();
        
        xray_image_files = all_xray_image_files{1};
        [landmark_tbls{dataset_indx}, joint_center_tbl] = LoadLandmarkFiles_foot(data_root_dir, x_ray_landmark_dir, patient_ID, view_direction, xray_image_files, x_ray_landmark_filename, bone_tbl.Row, leading_view, CT_file_name, CT_landmark_filename);
    end

    %%
    figure('Position',[300 300 [1600 900]*6/10], 'PaperPositionMode', 'auto', 'Color', 'w');
    for i=1:height(bone_tbl)
        subplot(2,2,i);
        auto_view1 = landmark_tbls{1}{1}.pos2D{i};
        auto_view2 = landmark_tbls{1}{2}.pos2D{i};
        manual_view1 = landmark_tbls{2}{1}.pos2D{i};
        manual_view2 = landmark_tbls{2}{2}.pos2D{i};
        landmark_error_2D_view1 = mean(sqrt((auto_view1(:,1,:)-manual_view1(:,1,:)).^2 + (auto_view1(:,2,:)-manual_view1(:,2,:)).^2),3);
        landmark_error_2D_view2 = mean(sqrt((auto_view2(:,1,:)-manual_view2(:,1,:)).^2 + (auto_view2(:,2,:)-manual_view2(:,2,:)).^2),3);
        error_3D_view1 = mean(sqrt(sum((landmark_tbls{1}{1}.pos3D{i}-landmark_tbls{2}{1}.pos3D{i}).^2,2)));
        error_3D_view2 = mean(sqrt(sum((landmark_tbls{1}{2}.pos3D{i}-landmark_tbls{2}{2}.pos3D{i}).^2,2)));
        plot(mean([landmark_error_2D_view1 landmark_error_2D_view2],2));
        set(gca,'xlim',[1 size(landmark_error_2D_view1,1)],'ylim',[0 15],'FontSize',14);
        xlabel('frame'); ylabel('error (pixel)');
        title(sprintf('%s, average error, 2D: %.3f pixel, 3D: %.3f mm',bone_tbl.Row{i}, ...
            mean([landmark_error_2D_view1(:); landmark_error_2D_view2(:)]),mean([error_3D_view1; error_3D_view2])),'FontSize',12);
    end
%     h=legend('x_trans','y_trans','z_trans','x_rot','y_rot','z_rot', 'interpreter','none','orientation','vertical','FontSize',10);
%     set(h,'Position',[0.95 0.6 0.01 0.3]);
    btitle(sprintf('ND%04d',patient_indx),18,[0 0 0], 'none', 0.98);
    print(gcf,'-dpng','-r300', fullfile(output_folder, sprintf('landmark_error_ND%04d_%s_%s.png',patient_indx, xray_image_files{1}, xray_image_files{2})));
end