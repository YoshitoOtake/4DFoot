addpath('C:\Users\otake\Documents\Programs\RegTools_build_VS2015\Matlab');
registration_results_folder = './registration_result_two_view_20201014_013955';

% registration_results_folder = 'registration_result_without_LCN_two_view_008';
root_dir = '../data';
patient_ID = 'ND0001';
Xp_dir = fullfile(root_dir, 'image/Xray', patient_ID);
CT_dir = fullfile(root_dir, 'image/CT', patient_ID);
Polygons_dir = fullfile(root_dir, 'polygon/CT', patient_ID);
Label_dir = fullfile(root_dir, 'image/CT', patient_ID);
load_bone_names_and_color_list; 
labels = {1, 3, 4, 5};
output_video_filename_prefix = '20201014_Foot2D3D';
GPU_IDs = [1]; %[0 1];

addpath('smoothpatch');
target_objects = [1 3 4 5];
bone_labels = {'tibia', 'talus', 'calc', 'navi'};
all_try_lambda = 0; %[1e-2 1e-3 1e-4 1e-5 1e-6 0];

fprintf('loading volumes of %s... ', patient_ID); timerID = tic;
[HU_img, HU_header] = mhdread( fullfile(CT_dir, 'image.mhd') );
[mask_img, mask_header] = mhdread( fullfile(CT_dir, 'label.mhd') );
volume_center_offset = eye(4);

myu_water = 0.2683;     % see RegTools.HU2Myu() function for detail
%         myu_img = RegTools.HU2Myu(HU_img, myu_water);
myu_img = HU2Myu_nonlinear(HU_img,myu_water); %max( (1000+single(HU_img)) * myu_water / 1000, 0 );
myu_img(myu_img<0) = 0;
mask_blur_width = 4; %3;
mask_img(1:floor(size(mask_img,1)/2),:,:) = 0;  % remove right foot

for erosion_size = 4 % [0 1 2 3 4 5 6 7 8 9 10]
    output_folder = fullfile(registration_results_folder, sprintf('erode_%03d_small_myu_water_test',erosion_size));
    if(~exist(output_folder, 'dir')), mkdir(output_folder); end
    for Lambda_L2_penalty = all_try_lambda(1)
    % for lambda_indx = 1:length(all_try_lambda)
    %     folder_names = {sprintf('Patient0000_%s_Lambda%0.0e', Fluoro_filename, all_try_lambda(lambda_indx))};
        view_direction = {'lateral', 'oblique'};
        xray_image_files = {'LL0027', 'LL0008'};
        folder_names = {fullfile(registration_results_folder, sprintf('%s_%s', patient_ID, xray_image_files{1}))};    % output_measurement_files = fullfile(Landmark_dir, strcat({folder_names.name}', '_measurement_results.csv'));

        % all_polygons_Xp_overlay = cell(length(folder_names), 1);
    %     volume_rendering_imgs = cell(length(folder_names), 2);
        Xp_imgs = cell(length(folder_names), 4);
        all_estimated_6Nxnum_frame = cell(length(folder_names), 1);
    %     all_estimated_wrt_global = cell(length(folder_names), 1);

        for folder_indx=1:length(folder_names)
            [img_2D_lateral, hdr_2D_lateral] = mhdread( fullfile(Xp_dir, view_direction{1}, sprintf('BP_%s.mhd',xray_image_files{1})) );
            [img_2D_oblique, hdr_2D_oblique] = mhdread( fullfile(Xp_dir, view_direction{2}, sprintf('BP_%s.mhd',xray_image_files{2})) );
            Xp_imgs(folder_indx,:) = {img_2D_lateral, hdr_2D_lateral, img_2D_oblique, hdr_2D_oblique};
            registration_result_file = fullfile( registration_results_folder, 'final_results.txt' );
            loaded_result = jsondecode(fileread(registration_result_file));
            [~, ~, num_bones, num_frames] = size(loaded_result.transformation_parameters);
            all_estimated_6Nxnum_frame{folder_indx} = reshape(RegTools.convert4x4ToTransRot_multi(reshape(loaded_result.transformation_parameters,4,4,[])),6,num_bones,num_frames);
        end

       %%
        camera_calib = jsondecode(fileread(fullfile(root_dir, 'calibration_files', sprintf('calibration_result_%s.json', patient_ID))));
        focal_length_pix = [camera_calib.camera1.focal_length camera_calib.camera2.focal_length]; %hdr_2D.SDD ./ hdr_2D.ElementSpacing(1:2);
        origin_2D_pix = [camera_calib.camera1.image_center camera_calib.camera2.image_center]; %[size(img_2D,1), size(img_2D,2)]./2;
        f1 = camera_calib.camera1.focal_length / mean(Xp_imgs{2}.ElementSpacing(1:2));
        f2 = camera_calib.camera2.focal_length / mean(Xp_imgs{4}.ElementSpacing(1:2));
        CameraPositions = cat(3, eye(4), camera_calib.relative_position');
        ProjectionMatrix1_pix = [[-f1, 0, camera_calib.camera1.image_center(1); 0, -f1, camera_calib.camera1.image_center(2); 0, 0, 1] [0 0 0]'];
        ProjectionMatrix2_pix = [[-f2, 0, camera_calib.camera2.image_center(1); 0, -f2, camera_calib.camera2.image_center(2); 0, 0, 1] [0 0 0]'] / CameraPositions(:,:,2);
        ProjectionMatrices_pix = cat(3, ProjectionMatrix1_pix, ProjectionMatrix2_pix);

        regTools = RegTools(GPU_IDs, [], 'log_file1.txt');
        all_DRRs = cell(length(folder_names), 1);
        all_DRRs_individual_bone = cell(length(folder_names), 1);
        for folder_indx=1:length(folder_names)
            if(isempty(all_estimated_6Nxnum_frame{folder_indx})), continue; end
            num_frames = size(all_estimated_6Nxnum_frame{folder_indx},3);
            [nu, nv, ~] = size(Xp_imgs{folder_indx,1});

            volumes = cell(length(labels), 1);
            for j=1:length(labels)
                mask = ismember(mask_img,labels{j});
                if(erosion_size>0),  mask = imerode(mask, strel('sphere',erosion_size)); end
%                 volumes{j} = myu_img.*mask;
%                 volumes{j} = ApplyEdgeBlurredMask(myu_img, mask, mask_blur_width);
                volumes{j} = ApplyErodedMask(myu_img, mask, mask_blur_width);
            end
            fprintf('done in %f sec\n', toc(timerID));

            fprintf('generating DRR of %s... ', patient_ID); timerID = tic;
%             num_frames = 10; all_estimated_6Nxnum_frame{folder_indx} = all_estimated_6Nxnum_frame{folder_indx}(:,:,1:num_frames);
            all_DRRs{folder_indx} = zeros([nu, nv, size(ProjectionMatrices_pix,3), num_frames], 'single');
            all_DRRs_individual_bone{folder_indx} = zeros([nu, nv, size(ProjectionMatrices_pix,3), num_frames, length(volumes)], 'single');
            planIDs = zeros(length(volumes),1,'int32');
            for i=1:length(volumes)
                planIDs(i) = regTools.CreateForwardProjectionPlan( volumes{i}, HU_header.ElementSpacing );
            end
            geomID = regTools.GenerateGeometry_3x4ProjectionMatrix( repmat(ProjectionMatrices_pix, [1 1 num_frames]), Xp_imgs{folder_indx,2}.ElementSpacing(1:2), [nu nv], [1 1] );
            regTools.SetLCNSigma(0);
            LCN_plan_ID = []; %regTools.CreateLCNComputationPlan([nu nv size(ProjectionMatrices_pix,3)], num_frames, 10);
            for i=1:length(volumes)
                transform = RegTools.convertTransRotTo4x4_multi(squeeze(all_estimated_6Nxnum_frame{folder_indx}(:,i,:)));
                transform = RegTools.MultiProd(transform, repmat(volume_center_offset,[1 1 num_frames])); % multiply one 4x4 matrix from the right
                all_DRRs_individual_bone{folder_indx}(:,:,:,:,i) = reshape(regTools.ForwardProject(planIDs(i), transform, [], 2, [], [], [], [], LCN_plan_ID),[nu, nv, size(ProjectionMatrices_pix,3), num_frames]);
                all_DRRs{folder_indx} = all_DRRs{folder_indx} + all_DRRs_individual_bone{folder_indx}(:,:,:,:,i);
            end
            for i=1:length(volumes)
                regTools.DeleteForwardProjectionPlan( planIDs(i) );
            end
            regTools.DeleteProjectionParametersArray( geomID );
            fprintf('done in %f sec\n', toc(timerID));
        end
        clear regTools;
        figure('Position',[100 100 [1200 600]*8/10], 'PaperPositionMode', 'auto', 'Color', 'w', 'InvertHardcopy', 'off');
%         clim = [-2 2];
        clim = [0 10];
        msubplot(1,2,1);
        im(reshape(permute(all_DRRs_individual_bone{folder_indx}(:,:,:,:,1),[1 2 4 3]),nu,nv,[])); colormap(gray(256)); colorbar; set(gca,'clim',clim);
        msubplot(1,2,2);
        im(reshape(permute(all_DRRs_individual_bone{folder_indx}(:,:,:,:,2),[1 2 4 3]),nu,nv,[])); colormap(gray(256)); colorbar; set(gca,'clim',clim);
        
        if(1)
        fprintf('saving images of %s... ', patient_ID); timerID = tic;
        for view=1:size(ProjectionMatrices_pix,3)
            for i=1:length(volumes)
                filename = fullfile(output_folder, [output_video_filename_prefix '_' sprintf('%s_%s_%s_Lambda%0.0e.mhd', patient_ID, xray_image_files{view}, bone_labels{i}, Lambda_L2_penalty)]);
                mhdwrite(squeeze(all_DRRs_individual_bone{folder_indx}(:,:,view,:,i)), filename, struct('ElementSpacing', Xp_imgs{2}.ElementSpacing, 'CompressedData', true));
            end
        end
        fprintf('done in %f sec\n', toc(timerID));
        end
        
        %%
        if(1)
    %         two_view_joined = cat(1, permute(squeeze(all_DRRs{1}(:,:,1,:)),[1 2 3]), permute(squeeze(all_DRRs{1}(:,:,2,:)),[1 2 3]));
    %         filename = fullfile(registration_results_folder, [output_video_filename_prefix '_' sprintf('%s_%s_Lambda%0.0e.mhd', PatientID, xray_image_files{1}, Lambda_L2_penalty)]);
    %         mhdwrite(two_view_joined, filename, struct('ElementSpacing', Xp_imgs{2}.ElementSpacing, 'CompressedData', true));
            four_objects = cell(4,2);
            for i=1:2
                for j=1:4
                    four_objects{j,i} = repmat(mat2gray(all_DRRs_individual_bone{1}(:,:,i,:,j)),[1 1 3 1]);
                end
            end
            
            four_objects_merged = cell(2,1);
            for i=1:2
                four_objects_merged{i} = repmat(mat2gray(all_DRRs{1}(:,:,i,:)),[1 1 3 1]);
                four_objects{1,i}(:,:,2:3,:) = 0;            four_objects{2,i}(:,:,[1 3],:) = 0;
                four_objects{3,i}(:,:,[1 2],:) = 0;          four_objects{4,i}(:,:,3,:) = 0;
            end
            four_object_joined = cell(2,1);
            for i=1:2
                four_object_joined{i} = cat(2, ...
                    cat(1, four_objects_merged{i} + four_objects{1,i}, four_objects_merged{i} + four_objects{2,i}), ...
                    cat(1, four_objects_merged{i} + four_objects{3,i}, four_objects_merged{i} + four_objects{4,i}) );
                four_object_joined{i}(nu+(-2:2),:,:,:) = 1.0; four_object_joined{i}(:,nv+(-2:2),:,:) = 1.0;
                four_object_joined{i}([1 2 end-1 end],:,:,:) = 1.0; four_object_joined{i}(:,[1 2 end-1 end],:,:) = 1.0;
            end
            filename1 = fullfile(registration_results_folder, [output_video_filename_prefix '_' sprintf('%s_%s_Lambda%0.0e.mhd', patient_ID, xray_image_files{1}, Lambda_L2_penalty)]);
            mhdwrite_NChannel(permute(four_object_joined{1},[1 2 4 3]), filename1, struct('ElementSpacing', Xp_imgs{2}.ElementSpacing, 'CompressedData', true));
            filename2 = fullfile(registration_results_folder, [output_video_filename_prefix '_' sprintf('%s_%s_Lambda%0.0e.mhd', patient_ID, xray_image_files{2}, Lambda_L2_penalty)]);
            mhdwrite_NChannel(permute(four_object_joined{2},[1 2 4 3]), filename2, struct('ElementSpacing', Xp_imgs{4}.ElementSpacing, 'CompressedData', true));
            return;
        end
    end
end