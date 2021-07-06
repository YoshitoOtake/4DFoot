output_folder = 'D:\Collaboration\Nara Medical University\20210218_DeepLabCut_analysis';
data_root_dir_landmark_auto = '\\scallop\user\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data_predicted'; 
data_root_dir_landmark_manual = '\\scallop\user\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data'; 

view_direction = {'lateral', 'oblique'};
x_ray_landmark_dir_auto = 'Xray_DLC_csv'; 
x_ray_landmark_dir_manual = 'Xray_csv'; 
x_ray_landmark_filename_auto = 'dlc_original.csv';
x_ray_landmark_filename_manual = 'CollectedData_miyamoto.csv'; 
all_landmarks = cell(5,2);
all_images = cell(5,2);
for patient_indx = 1:5
    patient_ID = sprintf('ND%04d', patient_indx);
    fprintf('loading all 2D landmarks of %s\n', patient_ID);
    for view=1:length(view_direction)
        all_files = dir( fullfile(data_root_dir_landmark_auto, 'landmark', x_ray_landmark_dir_auto, patient_ID, view_direction{view}, 'BP_*' ) );
        one_view_landmarks = cell(length(all_files),4);
        for j=1:length(all_files)
            landmark_2D_filename_auto = fullfile(all_files(j).folder, all_files(j).name, x_ray_landmark_filename_auto);
            [pos2D_auto, ~, likelihood] = LoadDeepLabCutCsvFile(landmark_2D_filename_auto);
            landmark_2D_filename_manual = fullfile(data_root_dir_landmark_manual, 'landmark', x_ray_landmark_dir_manual, patient_ID, view_direction{view}, all_files(j).name, x_ray_landmark_filename_manual);
            [pos2D_manual, ~] = LoadDeepLabCutCsvFile(landmark_2D_filename_manual);
            one_view_landmarks(j,:) = {pos2D_manual, pos2D_auto, likelihood', squeeze(sqrt(sum((pos2D_manual-pos2D_auto).^2,2)))'};
            
            if(j==1)
                img2D_fliename = fullfile(data_root_dir_landmark_manual, 'image', 'Xray', patient_ID, view_direction{view}, [all_files(j).name '.mhd']);
                [img2D, hdr2D] = mhdread( img2D_fliename );
            end
        end
        all_landmarks{patient_indx, view} = one_view_landmarks;
        show_frame = floor(size(img2D,3)/2);
        all_images{patient_indx, view} = {img2D(:,:,show_frame), show_frame};
    end
end

all_likelihood = cellfun(@(x) cell2mat(x(:,3)), all_landmarks, 'UniformOutput', false);
all_errors = cellfun(@(x) cell2mat(x(:,4)), all_landmarks, 'UniformOutput', false);

%%
figure('Position',[10 10 [1300 1000]*10/10], 'PaperPositionMode', 'auto', 'Color', 'w');
[hs,vs,tb,bb,lb,rb] = deal(0.05,0.05,0.05,0.05,0.05,0.05);
for patient_indx = 1:5
    patient_ID = sprintf('ND%04d', patient_indx);
    clf;
    for view=1:length(view_direction)
        for landmark_ID=1:size(all_likelihood{patient_indx,view},2)
            msubplot(4,6,(view-1)*12+landmark_ID,hs,vs,tb,bb,lb,rb);
            scatter( all_likelihood{patient_indx,view}(:,landmark_ID),  all_errors{patient_indx,view}(:,landmark_ID), '+' );
            xlabel('likelihood'); ylabel('Error (pixel)');
            set(gca, 'FontSize', 10, 'xlim', [0 1]);
            title(sprintf('view: %s, landmark %d', view_direction{view}, landmark_ID), 'FontSize', 12);
        end
    end
    btitle(sprintf('%s', patient_ID), 18, [0 0 0], 'none', 0.99);
    print(gcf, '-dpng', '-r200', fullfile(output_folder, sprintf('%s_DeepLabCut_likelihood_error.png', patient_ID)));
end

%%
figure('Position',[10 10 [600 1200]*10/10], 'PaperPositionMode', 'auto', 'Color', 'w'); colormap(gray(256));
[hs,vs,tb,bb,lb,rb] = deal(0.05,0.05,0.05,0.05,0.05,0.005);
for patient_indx = 1:5
    patient_ID = sprintf('ND%04d', patient_indx);
    clf;
    for view=1:2
        msubplot(2,1,view,hs,vs,tb,bb,lb,rb);
        imagesc(all_images{patient_indx,view}{1}'); axis image; axis off;
        hold on;
        show_frame = all_images{patient_indx,view}{2};
        landmark_coords = all_landmarks{patient_indx,view}{1,1}(:,:,show_frame);
        scatter(landmark_coords(:,1),landmark_coords(:,2),150,'r', '+','LineWidth',3);
        for j=1:size(landmark_coords,1)
            text(landmark_coords(j,1), landmark_coords(j,2), sprintf('%d',j), 'Color', 'y', 'FontSize', 24);
        end
        title(sprintf('%s', view_direction{view}), 'FontSize', 18);
    end
    btitle(sprintf('%s', patient_ID), 24, [0 0 0], 'none', 0.99);
    print(gcf, '-dpng', '-r200', fullfile(output_folder, sprintf('%s_DeepLabCut_likelihood_error_reference_images.png', patient_ID)));
end
