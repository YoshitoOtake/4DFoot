addpath(getenv('RegToolsPath'));
% auto_dir = '\\scallop\User\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data_predicted\registration_results\20210122_ver1';
% manual_dir = '\\scallop\User\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data\registration_results\20210122_ver1';
auto_dir = 'D:\Collaboration\Nara Medical University\20210214_ver1_Foot2D3D_registration_results_auto\merged';
manual_dir = 'D:\Collaboration\Nara Medical University\20210214_ver1_Foot2D3D_registration_results_manual\merged';
output_folder = auto_dir;

patient_indx_array = 4; % 1:3 ;%[1 2 3 4 5];
all_error_6xN_cell = cell(length(patient_indx_array),1);
fig1 = figure('Position',[300 300 [1600 800]*8/10], 'PaperPositionMode', 'auto', 'Color', 'w');
for patient_indx = patient_indx_array
    [patient_ID, start_frame, end_frame, all_xray_image_files, calibration_patient_ID] = patient_specific_setup(patient_indx);
    local_coordinate_file = fullfile('\\scallop\User\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data\local_coordinate\', patient_ID, 'local_coordinate.json');
    local_coordinate = jsondecode( fileread(local_coordinate_file) );
    for j=1:size(local_coordinate.center2local_4x4xN,3), local_coordinate.center2local_4x4xN(1:3,4,j) = local_coordinate.center2local_4x4xN(1:3,4,j) * -1; end    % wrt local coordinate
    
    auto_json = dir( fullfile(auto_dir, sprintf('ND%04d*.json',patient_indx)) );
    manual_json = dir( fullfile(manual_dir,sprintf('ND%04d*.json',patient_indx)) );

    auto_data = jsondecode( fileread(fullfile(auto_json(1).folder, auto_json(1).name)) );
    manual_data = jsondecode( fileread(fullfile(manual_json(1).folder, manual_json(1).name)) );

    all_error_trans = zeros( size(manual_data.transformation_parameters) );
    for i=1:size(manual_data.transformation_parameters,4)
        for j=1:size(manual_data.transformation_parameters,3)
            if(j<=size(local_coordinate.center2local_4x4xN,3)), local_coord = local_coordinate.center2local_4x4xN(:,:,j); else, local_coord = eye(4); end
            all_error_trans(:,:,j,i) = inv(manual_data.transformation_parameters(:,:,j,i)*local_coord)*auto_data.transformation_parameters(:,:,j,i)*local_coord;
%             all_error_trans(:,:,j,i) = manual_data.transformation_parameters(:,:,j,i);
        end
    end

    bone_list = jsondecode(fileread('ankle_new.json'));
    bone_names = cellfun(@(x) x{2},bone_list.color,'UniformOutput',false);
    bone_label = bone_names([1 3:14]); %{'tibia', 'talus', 'calcaneal', 'navicular', 'Medial_cuneiform', 'Intermediate_cuneiform', 'Lateral_cuneiform', 'cuboid', '1st_metatarsal', '}';

    all_error_6xN = zeros(6, size(all_error_trans,4), length(bone_label));
    for i=1:length(bone_label)
        all_error_6xN(:,:,i) = RegTools.convert4x4ToTransRot_multi(all_error_trans(:,:,i,:));
    end
    all_error_6xN_cell{patient_indx} = all_error_6xN;

    [hs,vs,tb,bb,lb,rb] = deal(0.05,0.07,0.14,0.05,0.05,0.05);
    for i=1:length(bone_label)
        msubplot(3,5,i, hs,vs,tb,bb,lb,rb);
        plot(manual_data.registration_frame, all_error_6xN_cell{patient_indx}(:,:,i)');
        set(gca,'xlim',[manual_data.registration_frame(1) manual_data.registration_frame(end)],'FontSize',10);
        set(gca,'ylim',[-5 5]);
        xlabel('frame'); ylabel('error (mm or deg)');
        mean_error = mean(abs(all_error_6xN_cell{patient_indx}(:,:,i)),2);
        title({sprintf('%s, MAE: %.3f, %.3f, %.3f', bone_label{i}, mean_error(1:3)), sprintf('%.3f, %.3f, %.3f',mean_error(4:6))},'FontSize',12,'interpreter','none');
    end
    h=legend('x_trans','y_trans','z_trans','x_rot','y_rot','z_rot', 'interpreter','none','orientation','horizontal','FontSize',10);
    set(h,'Position',[0.7 0.1 0.15 0.1]);
    btitle(sprintf('ND%04d',patient_indx),18,[0 0 0], 'none', 0.98);
    print(gcf,'-dpng','-r300', fullfile(output_folder, sprintf('evaluation_result_auto_vs_manual_ND%04d_all_bones.png',patient_indx)));
    
    axis_titles = {'X translation (mm)','Y translation (mm)','Z translation (mm)','X rotation (deg)','Y rotation (deg)','Z rotation (deg)',};
    num_frames = size(all_error_trans,4);
    for i=1:4 %length(bone_label)
        subplot(2,4,i);
        plot(manual_data.registration_frame, all_error_6xN_cell{patient_indx}(:,:,i)', 'LineWidth', 2);
        set(gca,'xlim',[manual_data.registration_frame(1) manual_data.registration_frame(end)],'FontSize',12);
        set(gca,'ylim',[-3 3]);
        xlabel('frame'); ylabel('error');
        title(sprintf('%s',bone_label{i}),'FontSize',20);
        if(i==1), h=legend(axis_titles, 'interpreter','none','orientation','horizontal','FontSize',12); set(h,'Position',[0.2 0.97 0.6 0.02]); end

        subplot(2,4,4+i);
        boxplot(abs(reshape(all_error_6xN_cell{patient_indx}(:,:,i),[],1)), repmat(axis_titles(:), num_frames, 1));
        set(gca,'ylim',[0 3],'FontSize',12,'XTickLabelRotation',90);
        ylabel('absolute error');
        title(sprintf('MAE: %.3f, %.3f, %.3f, %.3f, %.3f, %.3f',mean(abs(all_error_6xN_cell{patient_indx}(:,:,i)),2)),'FontSize',12);
    end
    print(gcf,'-dpng','-r300', fullfile(output_folder, sprintf('evaluation_result_auto_vs_manual_ND%04d_proximal_tarsal_bones.png',patient_indx)));
end

success_frames = cellfun(@(x) all(abs(x(:,:,1:4))<3,1),all_error_6xN_cell,'UniformOutput',false);
success_rates = cell2mat(cellfun(@(x) squeeze(sum(x,2)/size(x,2)*100), success_frames,'UniformOutput',false)');
for i=1:length(patient_indx_array)
    fprintf('ND%04d: %f, %f, %f, %f\n', patient_indx_array(i), success_rates(:,i));
end

