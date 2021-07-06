addpath(getenv('RegToolsPath'));
addpath('../util/smoothpatch');
addpath('../util');
% registration_results_root = 'D:\Collaboration\Nara Medical University\20210213_ver4_Foot2D3D_registration_results_auto\merged';
dataset_ID_array = {'auto', 'manual', 'semi-auto', 'bone_model_manual', 'bone_model_RSA'};
for dataset_indx = [1 2] %4:5
    dataset_ID = dataset_ID_array{dataset_indx} ;
    registration_results_root = fullfile('D:\Collaboration\Nara Medical University', sprintf('20210214_ver1_Foot2D3D_registration_results_%s', dataset_ID), 'merged');
    for patient_indx = 1:3 % -1 % [1 2 3 4 5]
        [patient_ID, start_frame, end_frame, all_xray_image_files, calibration_patient_ID] = patient_specific_setup(patient_indx);
        setup_dataset;

        registration_result_file = dir( fullfile(registration_results_root, sprintf('%s_*.json',patient_ID)) );
        loaded_result = jsondecode(fileread( fullfile(registration_result_file(1).folder, registration_result_file(1).name)));
         [~, output_video_filename_prefix, ~] = fileparts(registration_result_file(1).name);
        output_rendering_folder = fullfile(registration_results_root, 'renderings');

        % load data
        [ProjectionMatrices_pix_3D, DimSize_2D, camera_calib] = LoadCameraCalibration_foot(data_root_dir_image, calibration_patient_ID);
        mask_filter_width = 0; %4; %3;
        [myu_img, HU_header, mask_img] = Load3DImage_foot(data_root_dir_image, patient_ID, CT_file_name);
        xray_image_files = all_xray_image_files{1};
        [landmark_tbls, joint_center_tbl] = LoadLandmarkFiles_foot(data_root_dir_landmark, x_ray_landmark_dir, patient_ID, view_direction, xray_image_files, x_ray_landmark_filename, loaded_result.bone_names, leading_view, CT_file_name, CT_landmark_filename);
        label_2D_erosion = -1; % 5; %15;
        disable_log_correction = false;
        border_mask = 5;
        if(dataset_indx==5), med_filt_size = []; else, med_filt_size = [5 5]; end
        [img_2D_original, similarity_measure_mask] = Load2DImages_foot(data_root_dir_image, patient_ID, Xray_suffix, view_direction, strcat(xray_image_files, x_ray_image_file_suffix), leading_view, label_2D_erosion, disable_log_correction, border_mask, med_filt_size);
        LCN_sigma = 0; %0;
        img_2D_original = ApplyLCN_image_wise(img_2D_original, LCN_sigma);
        img_2D = img_2D_original(:,:,:,loaded_result.registration_frame);

        % load bone color
        bone_list = jsondecode(fileread('ankle_new.json'));
        first_bone_indx = cellfun(@(x) x(1), loaded_result.bone_label_indx);
        bone_colors = cell2mat( cellfun(@(x) [x{3} x{4} x{5}],bone_list.color(first_bone_indx),'UniformOutput',false) );

        % generate DRRs
        regTools = RegTools(0, [], 'log_file1.txt');
        DS = [1 1];
        [all_masked_volumes, all_ElementSpacing_volume, volume_offset_4x4xN] = prepareBoneVolumes(myu_img, mask_img, loaded_result.bone_label_indx, mask_filter_width, HU_header.ElementSpacing, regTools, DS);
        planIDs = zeros(length(all_masked_volumes),1,'int32');
        for i=1:length(all_masked_volumes)
            planIDs(i) = regTools.CreateForwardProjectionPlan( all_masked_volumes{i}, HU_header.ElementSpacing );
        end
        geomID = regTools.GenerateGeometry_3x4ProjectionMatrix( ProjectionMatrices_pix_3D, [1 1], DimSize_2D(1,:), [1 1] );

        all_DRRs_individual_bone = zeros([DimSize_2D(1,:), size(ProjectionMatrices_pix_3D,3), loaded_result.num_frames, length(all_masked_volumes)], 'single');
        fprintf('rendering DRR... frame');
        for k=1:loaded_result.num_frames
            fprintf(' %d',k);
            for i=1:length(all_masked_volumes)
                all_DRRs_individual_bone(:,:,:,k,i) = regTools.ForwardProject(planIDs(i), loaded_result.transformation_parameters(:,:,i,k)*volume_offset_4x4xN(:,:,i), [], 2);
            end
        end
        fprintf('\n');

        for i=1:length(all_masked_volumes)
            regTools.DeleteForwardProjectionPlan( planIDs(i) );
        end
        regTools.DeleteProjectionParametersArray( geomID );
        clear regTools;
        surface_patch = prepareSurfacePatches(all_masked_volumes(1,:), all_ElementSpacing_volume(1,:));

        %%
        figure('Position',[100 100 [1600 900]*8/10], 'PaperPositionMode', 'auto', 'Color', 'w', 'InvertHardcopy', 'off');
        colormap(flip(jet(256),1));
        [hs,vs,tb,bb,lb,rb] = deal(0.001, 0.001, 0.1, 0.01, 0.01, 0.01);

        % setting for side and top view 3D surface
        [hs2,vs2,tb2,bb2,lb2,rb2] = deal([0 0]+hs, [0 0]+vs, [0 0.1]+tb, [0 0]+bb, [0 0.20]+lb, [0 0]+rb);
        camera_up = [-1 0 0; 0 0 1];
        camera_pos = [0 0 -1; -1 0 0];
        camera_distance = [350 300];
            
        view_angles = [0 180 180 -90 90 270];
        num_column = 2;
        [sph_X, sph_Y, sph_Z] = sphere(8);
        sph_r = 6;
        sph_X = sph_X * sph_r; sph_Y = sph_Y * sph_r; sph_Z = sph_Z * sph_r;
        Xp_overlay_labels = [1 2]; % without patella
        % Xp_overlay_labels = [1 2 3 4 5 6 7 8 9]; % with patella
        DRR_clim = [0 20];
        TitleFontSize = 16;
        bottomTitleFontSize = 12;
        titleStr = { 'original Xray', 'Xray ROI', 'DRR', 'DRR ROI', 'overlay', 'overlay ROI' };
        titleYPos = 0.88;
        bottomYPos = 0.005;
        ROI_size = 100;

        volume_size_mm = HU_header.ElementSpacing.*HU_header.DimSize;
        origin_offset = RegTools.matTranslation(volume_size_mm/2.*[-1 -1 1]);

        video = VideoWriter(fullfile(registration_results_root, output_video_filename_prefix), 'MPEG-4');
        if(~exist(output_rendering_folder, 'dir')), mkdir(output_rendering_folder); end
        video.FrameRate = 10;
        open(video);

        % prepare grid for texture mapping
        XYZ_Source = cell(size(ProjectionMatrices_pix_3D,3),1);
        camera_array = {camera_calib.camera1, camera_calib.camera2};
        CameraPositions = cat(3, diag([-1 1 1 1]), diag([-1 1 1 1])*camera_calib.relative_position');
        for view=1:size(ProjectionMatrices_pix_3D,3)
            DimSize = camera_array{view}.dim_size;
            Xp_half_size = DimSize.*camera_array{view}.element_spacing([1 2])/2;
            ImageCenter_mm = camera_array{view}.image_center([1 2]).*camera_array{view}.element_spacing([1 2]);
    %             [X,Y] = meshgrid((0:(DimSize(1)-1))*camera_array{view}.element_spacing(1),((DimSize(2)-1):-1:0)*camera_array{view}.element_spacing(2));
            [X,Y] = meshgrid((0:(DimSize(1)-1))*camera_array{view}.element_spacing(1),(0:(DimSize(2)-1))*camera_array{view}.element_spacing(2));
            Z = zeros(size(X));
            XYZ_Source{view} = CameraPositions(:,:,view) * [eye(3) [-ImageCenter_mm; -camera_array{view}.focal_length]; 0 0 0 1] * [X(:)'; Y(:)'; Z(:)'; ones(1,numel(X))];
        end

        [nu, nv] = deal(DimSize_2D(1,1),DimSize_2D(1,2));
        base_bone_indx = 2;
        base_bone_transform = zeros(4,4,loaded_result.num_frames);
        for k=1:loaded_result.num_frames
            base_bone_transform(:,:,k) = loaded_result.transformation_parameters(:,:,base_bone_indx,k)*volume_offset_4x4xN(:,:,base_bone_indx);
        end
        base_bone_transform = [eye(3) mean(base_bone_transform(1:3,4,:),3); 0 0 0 1];

        for k=1:loaded_result.num_frames
            % calculate default projection matrix
            all_2Dlandmark_center = [   
                                        median(cell2mat(cellfun(@(x) permute(x(loaded_result.registration_frame(k),:,:),[3 2 1]), landmark_tbls{1,1}.pos2D, 'UniformOutput', false))), ...
                                        median(cell2mat(cellfun(@(x) permute(x(loaded_result.registration_frame(k),:,:),[3 2 1]), landmark_tbls{2,1}.pos2D, 'UniformOutput', false)))];
            all_2Dlandmark_center = min(max(all_2Dlandmark_center,ROI_size+1),[nu-ROI_size-1, nv-ROI_size-1, nu-ROI_size-1, nv-ROI_size-1]);

            ROI_corners = zeros(2,4);
            ROI_corners(1,:) = floor(all_2Dlandmark_center([1 1 2 2]) + ROI_size .* [-1 1 -1 1]); 
            ROI_corners(2,:) = floor(all_2Dlandmark_center([3 3 4 4]) + ROI_size .* [-1 1 -1 1]); 

            clf;

            for v=1:2
                if(v==1), render_surface = 1:length(surface_patch); else, render_surface = 2:length(surface_patch); end
                msubplot(3,5,12+v,hs2(v),vs2(v),tb2(v),bb2(v),lb2(v),rb2(v));
                hold on;
                for i=render_surface
                    vert = ApplyTransform4x4(surface_patch{i}.vertices, inv(base_bone_transform)*loaded_result.transformation_parameters(:,:,i,k)*volume_offset_4x4xN(:,:,i)); %diag([1 -1 1 1])*RegTools.convertTransRotTo4x4(all_estimated_6Nxnum_frame{folder_indx}(:,i,k)));
                    patch('Faces', surface_patch{i}.faces, 'Vertices', vert, 'FaceColor', bone_colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 1.0); %object_alphas(indx_in_polygon)); %, 'EdgeColor', 'None', 'CData', u_disc, 'FaceColor', 'interp');
                end
                set(gca, 'clim', [0 10]);
                camera_target = [0 0 0]; %base_bone_transform(1:3,4)'; %all_moved3Dlandmark_center(1,:) + [200 0 0];
                camproj('perspective'); camtarget(camera_target); camup(camera_up(v,:)); axis off equal;
                set(gca,'CameraViewAngleMode', 'manual', 'CameraViewAngle', 20);
                material dull;
                campos(camera_target+camera_pos(v,:)*camera_distance(v)); 
                camlight('headlight');
            end
            
            for view=1:size(ProjectionMatrices_pix_3D,3)
                fixed = img_2D(:,:,view,k);
                fixed = mat2gray(fixed, fixed_clim); %double([min(fixed(:)) max(fixed(:))]));
                floating = sum(all_DRRs_individual_bone(:,:,view,k,:),5);

                msubplot(3,6,(view-1)*6+1,hs,vs,tb,bb,lb,rb);
                imc(repmat(fixed',[1 1 3])); axis off;
                line(ROI_corners(view,[1 2 2 1 1]), ROI_corners(view,[3 3 4 4 3]), 'Color', 'y', 'LineWidth', 1);
                if(view==1), title(titleStr{1}, 'FontSize', TitleFontSize); end

                msubplot(3,6,(view-1)*6+2,hs,vs,tb,bb,lb,rb);
                imc(repmat(mat2gray(fixed(ROI_corners(view,1):ROI_corners(view,2),ROI_corners(view,3):ROI_corners(view,4))'),[1 1 3])); axis off;
                line([0 0 1 1 0]*ROI_size*2, [0 1 1 0 0]*ROI_size*2, 'Color', 'y', 'LineWidth', 2);
                if(view==1), title(titleStr{2}, 'FontSize', TitleFontSize); end

                msubplot(3,6,(view-1)*6+3,hs,vs,tb,bb,lb,rb);
                imc(repmat(mat2gray(floating',moving_clim),[1 1 3])); axis off;
                line(ROI_corners(view,[1 2 2 1 1]), ROI_corners(view,[3 3 4 4 3]), 'Color', 'y', 'LineWidth', 1);
                if(view==1), title(titleStr{3}, 'FontSize', TitleFontSize); end

                msubplot(3,6,(view-1)*6+4,hs,vs,tb,bb,lb,rb);
                imc(repmat(mat2gray(floating(ROI_corners(view,1):ROI_corners(view,2),ROI_corners(view,3):ROI_corners(view,4))',moving_clim),[1 1 3])); axis off;
                line([0 0 1 1 0]*ROI_size*2, [0 1 1 0 0]*ROI_size*2, 'Color', 'y', 'LineWidth', 2);
                if(view==1), title(titleStr{4}, 'FontSize', TitleFontSize); end

                msubplot(3,6,(view-1)*6+5,hs,vs,tb,bb,lb,rb);
                fixed_rgb = permute( repmat( fixed, [1 1 1 3]), [1 2 4 3]);
                contour = 255 - canny('image',uint8(floating./max(floating(:)).*255), 'thigh', 0.5);
                contour = imdilate(contour, strel('square', 2));
                fixed_rgb(:,:,1) = fixed_rgb(:,:,1) + double(logical(contour));
                imc( permute(fixed_rgb,[2 1 3 4])); axis off;
                line(ROI_corners(view,[1 2 2 1 1]), ROI_corners(view,[3 3 4 4 3]), 'Color', 'y', 'LineWidth', 1);
                if(view==1), title(titleStr{5}, 'FontSize', TitleFontSize); end

                msubplot(3,6,(view-1)*6+6,hs,vs,tb,bb,lb,rb);
                imc( permute(fixed_rgb(ROI_corners(view,1):ROI_corners(view,2),ROI_corners(view,3):ROI_corners(view,4),:,:),[2 1 3 4])); axis off;
                line([0 0 1 1 0]*ROI_size*2, [0 1 1 0 0]*ROI_size*2, 'Color', 'y', 'LineWidth', 2);
                if(view==1), title(titleStr{6}, 'FontSize', TitleFontSize); end

                msubplot(3,5,10+view,hs,vs,tb,bb,lb,rb);
                hold on;
                surface(reshape(XYZ_Source{view}(1,:),size(X)),reshape(XYZ_Source{view}(2,:),size(X)),reshape(XYZ_Source{view}(3,:),size(X)), ...
                    'FaceColor','texturemap','EdgeColor','none', 'CData',repmat(mat2gray(img_2D(:,:,view,k)',fixed_clim),[1 1 3]), 'CDataMapping', 'direct', 'FaceLighting', 'none');
                for i=1:length(surface_patch)
                    vert = ApplyTransform4x4(surface_patch{i}.vertices, CameraPositions(:,:,1)*loaded_result.transformation_parameters(:,:,i,k)*volume_offset_4x4xN(:,:,i)); %diag([1 -1 1 1])*RegTools.convertTransRotTo4x4(all_estimated_6Nxnum_frame{folder_indx}(:,i,k)));
                    patch('Faces', surface_patch{i}.faces, 'Vertices', vert, 'FaceColor', bone_colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5); %object_alphas(indx_in_polygon)); %, 'EdgeColor', 'None', 'CData', u_disc, 'FaceColor', 'interp');
                end
                camproj('perspective');
                campos(CameraPositions(1:3,4,view)); 
                camtarget(ApplyTransform4x4([0 0 -camera_array{view}.focal_length],CameraPositions(:,:,view))); 
                camup(-CameraPositions(1:3,2,view)); 
                axis off equal;
                set(gca,'CameraViewAngleMode', 'manual', 'CameraViewAngle', rad2deg(atan2(min(Xp_half_size),camera_array{view}.focal_length))*2);
                camlight('headlight');
            end

            btitle(sprintf('%s, %s & %s, Frame %d', patient_ID, xray_image_files{:}, loaded_result.registration_frame(k)), 24, [0 0 0], 'tex', 0.97);
            drawnow;
            frame = getframe(gcf);
            writeVideo(video, frame);
        end
        close(video);
    end   
end