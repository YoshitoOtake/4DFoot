addpath(getenv('RegToolsPath'));
% predict_dir = '\\scallop\user\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data_model\registration_results\20210203_MODEL1_ver4';
% gt_dir = '\\scallop\User\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data_model\registration_results\20210207_RSA_ver3';
predict_dir = 'D:\Collaboration\Nara Medical University\Foot2D3D_registration_results\20210213_ver4_Foot2D3D_registration_results_bone_model_manual\merged';
gt_dir = 'D:\Collaboration\Nara Medical University\Foot2D3D_registration_results\20210213_ver4_Foot2D3D_registration_results_bone_model_RSA\merged';
local_coordinate_file = '\\scallop\User\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data_model\local_coordinate\MODEL1\local_coordinate.json';
local_coordinate = jsondecode( fileread(local_coordinate_file) );
for j=1:size(local_coordinate.center2local_4x4xN,3), local_coordinate.center2local_4x4xN(1:3,4,j) = local_coordinate.center2local_4x4xN(1:3,4,j) * -1; end    % wrt local coordinate
% local_coordinate.center2local_4x4xN = repmat(eye(4), [1 1 size(local_coordinate.center2local_4x4xN,3)]);    % wrt global coordinate

predict_json = dir( fullfile(predict_dir, '*.json') );
gt_json = dir( fullfile(gt_dir, '*.json') );

predict_data = jsondecode( fileread(fullfile(predict_json(1).folder, predict_json(1).name)) );
predict_data.volume_offset_4x4xN = reshape(predict_data.volume_offset_4x4xN, 4, 4, predict_data.num_transforms);
gt_data = jsondecode( fileread(fullfile(gt_json(1).folder, gt_json(1).name)) );

output_folder = predict_dir;
all_error_trans = zeros( size(gt_data.transformation_parameters) );
for i=1:size(gt_data.transformation_parameters,4)
    for j=1:size(gt_data.transformation_parameters,3)
        all_error_trans(:,:,j,i) = inv(gt_data.transformation_parameters(:,:,j,i)*local_coordinate.center2local_4x4xN(:,:,j))*predict_data.transformation_parameters(:,:,j,i)*local_coordinate.center2local_4x4xN(:,:,j);
%         all_error_trans(:,:,j,i) = predict_data.transformation_parameters(:,:,j,i);
%         all_error_trans(:,:,j,i) = gt_data.transformation_parameters(:,:,j,i);
    end
end

bone_list = jsondecode(fileread('ankle_new.json'));
bone_names = cellfun(@(x) x{2},bone_list.color,'UniformOutput',false);
bone_label = bone_names([1 3:5]); %{'tibia', 'talus', 'calcaneal', 'navicular', 'Medial_cuneiform', 'Intermediate_cuneiform', 'Lateral_cuneiform', 'cuboid', '1st_metatarsal', '}';

%%
figure('Position',[300 300 [1600 900]*10/10], 'PaperPositionMode', 'auto', 'Color', 'w');
axis_titles = {'X translation (mm)','Y translation (mm)','Z translation (mm)','X rotation (deg)','Y rotation (deg)','Z rotation (deg)',};
num_frames = size(all_error_trans,4);
all_error_6xN = zeros(6, num_frames, length(bone_label));
for i=1:length(bone_label)
    subplot(2,4,i);
    all_error_6xN(:,:,i) = RegTools.convert4x4ToTransRot_multi(all_error_trans(:,:,i,:));
    plot(all_error_6xN(:,:,i)', 'LineWidth', 2);
    set(gca,'xlim',[1 size(all_error_6xN(:,:,i),2)],'FontSize',12);
    set(gca,'ylim',[-3 3]);
    xlabel('frame'); ylabel('error');
    title(sprintf('%s',bone_label{i}),'FontSize',20);
    if(i==1), h=legend(axis_titles, 'interpreter','none','orientation','horizontal','FontSize',12); set(h,'Position',[0.2 0.97 0.6 0.02]); end

    subplot(2,4,4+i);
    boxplot(abs(reshape(all_error_6xN(:,:,i),[],1)), repmat(axis_titles(:), num_frames, 1));
    set(gca,'ylim',[0 3],'FontSize',12,'XTickLabelRotation',90);
    ylabel('absolute error');
    title(sprintf('MAE: %.3f, %.3f, %.3f, %.3f, %.3f, %.3f',mean(abs(all_error_6xN(:,:,i)),2)),'FontSize',12);
end
fprintf('absolute translation error of tarsal bones: %f +- %f, absolute rotation error of tarsal bones: %f +- %f\n', ...
    mean(reshape(abs(all_error_6xN(1:3,:,2:4)),[],1)), std(reshape(abs(all_error_6xN(1:3,:,2:4)),[],1)), ...
    mean(reshape(abs(all_error_6xN(4:6,:,2:4)),[],1)), std(reshape(abs(all_error_6xN(4:6,:,2:4)),[],1)) );
print(gcf,'-dpng','-r300', fullfile(output_folder, 'evaluation_result_with_bone_model.png'));