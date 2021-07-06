addpath('C:\Users\otake\Documents\Programs\RegTools_build_VS2015\Matlab');
% registration_results_folder = 'D:\Collaboration\Nara Medical University\20200910_Foot2D3D\Matlab_script/results_ND0001_LL0027_LL0008_20201027_065914_success';

% root_dir = '../data';
root_dir = 'D:\Collaboration\Nara Medical University\20200910_Foot2D3D\data';
% all_try_lambda = [1e0 1e-1 1e-2 1e-3 1e-4 1e-5 0];
patient_ID = 'ND0001';
Xp_dir = fullfile(root_dir, 'image/Xray', patient_ID);
CT_dir = fullfile(root_dir, 'image/CT', patient_ID);
Label_dir = fullfile(root_dir, 'image/CT', patient_ID);
Polygons_dir = fullfile(root_dir, 'polygon');
Landmark_dir = fullfile(root_dir, 'landmark/CT', patient_ID);
load_bone_names_and_color_list; 
labels = {[1 2], 3, 4, 5, 6, 7, 8, 9}; %, 10, 11, 12, 13, 14};
output_video_filename_prefix = '20201016_Foot2D3D';
output_rendering_folder = '20201016_polygon_visualization';
all_xray_image_files = { {'LL0027', 'LL0008'}, {'LL0028', 'LL0009'}, {'LL0029', 'LL0010'}, {'LL0034', 'LL0016'}, {'LL0035', 'LL0017'}, {'LL0036', 'LL0018'}, {'LL0037', 'LL0019'}};

addpath('../util/smoothpatch');
addpath('../util');
target_objects = [1 3:9]; %[1 3:14];
all_try_lambda = [0]; %[1e-2 1e-3 1e-4 1e-5 1e-6 0];
% for Lambda_L2_penalty = all_try_lambda(1)
for experiment_ID = 1:7
% for lambda_indx = 1:length(all_try_lambda)
%     folder_names = {sprintf('Patient0000_%s_Lambda%0.0e', Fluoro_filename, all_try_lambda(lambda_indx))};
    view_direction = {'lateral', 'oblique'};
    xray_image_files = all_xray_image_files{experiment_ID};
    leading_view = 2;   
    registration_results_root = 'D:\Collaboration\Nara Medical University\20200910_Foot2D3D\Matlab_script/';
    registration_results_folder = fullfile( registration_results_root, sprintf('results_ND0001_%s_%s_20201027_221600_success', xray_image_files{:}) );
%     registration_results_folder = 'D:\Collaboration\Nara Medical University\20200910_Foot2D3D\Matlab_script/results_ND0001_LL0027_LL0008_20201027_065914_success';
    folder_names = {patient_ID}; %{fullfile(registration_results_root, registration_results_folder, sprintf('%s_%s', patient_ID, xray_image_files{1}))};    % output_measurement_files = fullfile(Landmark_dir, strcat({folder_names.name}', '_measurement_results.csv'));

    volume_rendering_imgs = cell(length(folder_names), 2);
    Xp_imgs = cell(length(folder_names), 2);
    Landmarks = cell(length(folder_names), 1);
    all_estimated_6Nxnum_frame = cell(length(folder_names), 1);
    all_estimated_wrt_global = cell(length(folder_names), 1);
    all_landmarks_wrt_global = cell(length(folder_names), 1);

    all_polygons = cell(length(target_objects),1);
    for j=1:length(target_objects)
        fprintf('%s of %s (%d/%d)...', object_names{target_objects(j)}, patient_ID, j, length(target_objects)); timerID_polygon = tic;
        input_ply_filename = fullfile(Polygons_dir, patient_ID, sprintf('%s-%s.ply',patient_ID, object_names{target_objects(j)}));
        [f, v] = ply_readEx(input_ply_filename,'tri');
        all_polygons{j} = struct('faces',int32(f'),'vertices',ApplyTransform4x4(v', diag([-1 -1 1 1])));  % Slicer coordinate -> Matlab coordinate
        fprintf('done in %f sec\n', toc(timerID_polygon));
    end
    
    for folder_indx=1:length(folder_names)
        label_2D_erosion = 0; %15;
        [img_2D_original, similarity_measure_mask] = Load2DImages_foot(root_dir, patient_ID, view_direction, xray_image_files, leading_view, label_2D_erosion, true);
        Xp_imgs(folder_indx,:) = {squeeze(img_2D_original(:,:,1,:)), squeeze(img_2D_original(:,:,2,:))};
        registration_result_file = fullfile( registration_results_folder, 'final_results.json' );
        loaded_result = jsondecode(fileread(registration_result_file));
        [~, ~, num_bones, num_frames] = size(loaded_result.transformation_parameters);
        all_estimated_6Nxnum_frame{folder_indx} = reshape(RegTools.convert4x4ToTransRot_multi(reshape(loaded_result.transformation_parameters,4,4,[])),6,num_bones,num_frames);
    end
    [ProjectionMatrices_pix_3D, ~, camera_calib] = LoadCameraCalibration_foot(root_dir, patient_ID);
    camera_calib.camera_array = {camera_calib.camera1, camera_calib.camera2};
    CameraPositions = cat(3, eye(4), camera_calib.relative_position');
    
   %%
    regTools = RegTools(0, [], 'log_file1.txt');
    all_DRRs = cell(length(folder_names), 1);
    all_DRRs_individual_bone = cell(length(folder_names), 1);
    for folder_indx=1:length(folder_names)
        if(isempty(all_estimated_6Nxnum_frame{folder_indx})), continue; end

        % calculate default projection matrix
        num_frames = size(all_estimated_6Nxnum_frame{folder_indx},3);
        [nu, nv, ~] = size(Xp_imgs{folder_indx,1});

        fprintf('loading volumes of %s... ', patient_ID); timerID = tic;
        mask_filter_width = 0; %4; %3;
        [myu_img, HU_header, mask_img] = Load3DImage_foot(root_dir, patient_ID);
        mask_blur_width = 4; %3;
        volumes = cell(length(labels), 1);
        for j=1:length(labels)
            mask = ismember(mask_img,labels{j});
            volumes{j} = ApplyEdgeBlurredMask(myu_img, mask, mask_filter_width);
%             volumes{j} = ApplyErodedMask(myu_img, mask, mask_filter_width);
%             volumes{j} = myu_img.*mask;
        end
        fprintf('done in %f sec\n', toc(timerID));

        fprintf('generating DRR of %s... ', patient_ID); timerID = tic;
        all_DRRs{folder_indx} = zeros([nu, nv, size(ProjectionMatrices_pix_3D,3), num_frames], 'single');
        all_DRRs_individual_bone{folder_indx} = zeros([nu, nv, size(ProjectionMatrices_pix_3D,3), num_frames, length(volumes)], 'single');
        planIDs = zeros(length(volumes),1,'int32');
        for i=1:length(volumes)
            planIDs(i) = regTools.CreateForwardProjectionPlan( volumes{i}, HU_header.ElementSpacing );
        end
        geomID = regTools.GenerateGeometry_3x4ProjectionMatrix( ProjectionMatrices_pix_3D, [1 1], [nu nv], [1 1] );
        for k=1:num_frames
            for i=1:length(volumes)
                transform = RegTools.convertTransRotTo4x4(all_estimated_6Nxnum_frame{folder_indx}(:,i,k));
                all_DRRs_individual_bone{folder_indx}(:,:,:,k,i) = regTools.ForwardProject(planIDs(i), transform, [], 2);
                all_DRRs{folder_indx}(:,:,:,k) = all_DRRs{folder_indx}(:,:,:,k) + all_DRRs_individual_bone{folder_indx}(:,:,:,k,i);
            end
        end
        for i=1:length(volumes)
            regTools.DeleteForwardProjectionPlan( planIDs(i) );
        end
        regTools.DeleteProjectionParametersArray( geomID );
        fprintf('done in %f sec\n', toc(timerID));
        
    end
    clear regTools;
    if(0)
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
        mhdwrite_NChannel(permute(four_object_joined{1},[1 2 4 3]), filename1, struct('ElementSpacing', camera_calib.camera1.element_spacing, 'CompressedData', true));
        filename2 = fullfile(registration_results_folder, [output_video_filename_prefix '_' sprintf('%s_%s_Lambda%0.0e.mhd', patient_ID, xray_image_files{2}, Lambda_L2_penalty)]);
        mhdwrite_NChannel(permute(four_object_joined{2},[1 2 4 3]), filename2, struct('ElementSpacing', camera_calib.camera2.element_spacing, 'CompressedData', true));
    end

    %%
    figure('Position',[100 100 [1600 600]*8/10], 'PaperPositionMode', 'auto', 'Color', 'w', 'InvertHardcopy', 'off');
    colormap(flip(jet(256),1));
    [hs,vs,tb,bb,lb,rb] = deal(0.001, 0.001, 0.1, 0.01, 0.01, 0.01);
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
    titleStr = { 'original fluoro', 'fluoro ROI', 'DRR', 'DRR ROI', 'overlay', 'overlay ROI' };
    titleYPos = 0.88;
    bottomYPos = 0.005;
    
    volume_size_mm = HU_header.ElementSpacing.*HU_header.DimSize;
    origin_offset = RegTools.matTranslation(volume_size_mm/2.*[-1 -1 1]);
    bone_label = object_names(target_objects)';
    landmark_tbls = LoadLandmarkFiles_foot(root_dir, patient_ID, view_direction, xray_image_files, bone_label, leading_view);
    
    for folder_indx=1:length(folder_names)
%         v = VideoWriter(fullfile(registration_results_folder, [output_video_filename_prefix '_' sprintf('%s_%s_Lambda%0.0e.mp4', patient_ID, xray_image_files{1}, Lambda_L2_penalty)]), 'MPEG-4');
        v = VideoWriter(fullfile(registration_results_folder, [output_video_filename_prefix '_' sprintf('%s_%s_%s_overlay_video.mp4', patient_ID, xray_image_files{:})]), 'MPEG-4');
        output_rendering_folder_full = fullfile(folder_names{folder_indx}, output_rendering_folder);
        if(~exist(output_rendering_folder_full, 'dir')), mkdir(output_rendering_folder_full); end
        v.FrameRate = 10;
        open(v);

        [nu, nv, ~] = size(Xp_imgs{folder_indx,1});
        num_frames = size(all_estimated_6Nxnum_frame{folder_indx},3);
        for k=1:num_frames
            % calculate default projection matrix
            all_2Dlandmark_center = [   mean(cell2mat(cellfun(@(x) permute(x(k,:,:),[3 2 1]), landmark_tbls{1,1}.pos2D, 'UniformOutput', false))), ...
                                        mean(cell2mat(cellfun(@(x) permute(x(k,:,:),[3 2 1]), landmark_tbls{2,1}.pos2D, 'UniformOutput', false)))];
            ROI_size = 100;
            all_2Dlandmark_center = min(max(all_2Dlandmark_center,ROI_size+1),[nu-ROI_size-1, nv-ROI_size-1, nu-ROI_size-1, nv-ROI_size-1]);
            
            ROI_corners = zeros(2,4);
            ROI_corners(1,:) = floor(all_2Dlandmark_center([1 1 2 2]) + ROI_size .* [-1 1 -1 1]); 
            ROI_corners(2,:) = floor(all_2Dlandmark_center([3 3 4 4]) + ROI_size .* [-1 1 -1 1]); 
            
            clf;
            for view=1:size(ProjectionMatrices_pix_3D,3)
                fixed = Xp_imgs{folder_indx,view}(:,:,k);
                fixed = mat2gray(fixed,double([min(fixed(:)) max(fixed(:))]));
                floating = all_DRRs{folder_indx}(:,:,view,k);
                
                msubplot(2,6,(view-1)*6+1,hs,vs,tb,bb,lb,rb);
                imc(repmat(fixed',[1 1 3])); axis off;
                line(ROI_corners(view,[1 2 2 1 1]), ROI_corners(view,[3 3 4 4 3]), 'Color', 'y', 'LineWidth', 1);
                if(view==1), title(titleStr{1}, 'FontSize', TitleFontSize); end

                msubplot(2,6,(view-1)*6+2,hs,vs,tb,bb,lb,rb);
                imc(repmat(mat2gray(fixed(ROI_corners(view,1):ROI_corners(view,2),ROI_corners(view,3):ROI_corners(view,4))'),[1 1 3])); axis off;
                line([0 0 1 1 0]*ROI_size*2, [0 1 1 0 0]*ROI_size*2, 'Color', 'y', 'LineWidth', 2);
                if(view==1), title(titleStr{2}, 'FontSize', TitleFontSize); end

                msubplot(2,6,(view-1)*6+3,hs,vs,tb,bb,lb,rb);
                imc(repmat(mat2gray(all_DRRs{folder_indx}(:,:,view,k)',DRR_clim),[1 1 3])); axis off;
                line(ROI_corners(view,[1 2 2 1 1]), ROI_corners(view,[3 3 4 4 3]), 'Color', 'y', 'LineWidth', 1);
                if(view==1), title(titleStr{3}, 'FontSize', TitleFontSize); end

                msubplot(2,6,(view-1)*6+4,hs,vs,tb,bb,lb,rb);
                imc(repmat(mat2gray(all_DRRs{folder_indx}(ROI_corners(view,1):ROI_corners(view,2),ROI_corners(view,3):ROI_corners(view,4),view,k)',DRR_clim),[1 1 3])); axis off;
                line([0 0 1 1 0]*ROI_size*2, [0 1 1 0 0]*ROI_size*2, 'Color', 'y', 'LineWidth', 2);
                if(view==1), title(titleStr{4}, 'FontSize', TitleFontSize); end

                msubplot(2,6,(view-1)*6+5,hs,vs,tb,bb,lb,rb);
                fixed_rgb = permute( repmat( fixed, [1 1 1 3]), [1 2 4 3]);
                contour = 255 - canny('image',uint8(floating./max(floating(:)).*255), 'thigh', 0.5);
                contour = imdilate(contour, strel('square', 2));
                fixed_rgb(:,:,1) = fixed_rgb(:,:,1) + double(logical(contour));
                imc( permute(fixed_rgb,[2 1 3 4])); axis off;
                line(ROI_corners(view,[1 2 2 1 1]), ROI_corners(view,[3 3 4 4 3]), 'Color', 'y', 'LineWidth', 1);
                if(view==1), title(titleStr{5}, 'FontSize', TitleFontSize); end

                msubplot(2,6,(view-1)*6+6,hs,vs,tb,bb,lb,rb);
                imc( permute(fixed_rgb(ROI_corners(view,1):ROI_corners(view,2),ROI_corners(view,3):ROI_corners(view,4),:,:),[2 1 3 4])); axis off;
                line([0 0 1 1 0]*ROI_size*2, [0 1 1 0 0]*ROI_size*2, 'Color', 'y', 'LineWidth', 2);
                if(view==1), title(titleStr{6}, 'FontSize', TitleFontSize); end
            end
            
            XYZ_Source = cell(size(ProjectionMatrices_pix_3D,3),1);
            
            for view=1:size(ProjectionMatrices_pix_3D,3)
                DimSize = camera_calib.camera_array{view}.dim_size;
                Xp_half_size = DimSize.*camera_calib.camera_array{view}.element_spacing([1 2])/2;
                [X,Y] = meshgrid((0:(DimSize(1)-1))*camera_calib.camera_array{view}.element_spacing(1),((DimSize(1)-1):-1:0)*camera_calib.camera_array{view}.element_spacing(2));
                Z = zeros(size(X));
                XYZ_Source{view} = inv(CameraPositions(:,:,view)) * [eye(3) [-Xp_half_size(:); -camera_calib.camera_array{view}.focal_length]; 0 0 0 1] * [X(:)'; Y(:)'; Z(:)'; ones(1,numel(X))];
                camtarget_Z = camera_calib.camera_array{view}.focal_length;
            end
if(0)
            msubplot(3,4,9,hs,vs,tb,bb,lb,rb);
%             msubplot(1,2,1);
            hold on;
            surface(reshape(XYZ_Source{1}(1,:),size(X)),reshape(XYZ_Source{1}(2,:),size(X)),reshape(XYZ_Source{1}(3,:),size(X)), ...
                'FaceColor','texturemap','EdgeColor','none', 'CData',repmat(mat2gray(Xp_imgs{folder_indx,1}(:,:,k)'),[1 1 3]), 'CDataMapping', 'direct', 'FaceLighting', 'none');
            for i=1:length(target_objects)
                vert = ApplyTransform4x4(all_polygons{i}.vertices, diag([1 -1 1 1])*RegTools.convertTransRotTo4x4(all_estimated_6Nxnum_frame{folder_indx}(:,i,k)));
                patch('Faces', all_polygons{i}.faces, 'Vertices', vert, 'FaceColor', object_colors(target_objects(i),:), 'EdgeColor', 'none', 'FaceAlpha', 0.3); %object_alphas(indx_in_polygon)); %, 'EdgeColor', 'None', 'CData', u_disc, 'FaceColor', 'interp');
            end
%             draw_coordinate_system(eye(4), 50, 3);
            camproj('perspective'); camtarget([0 0 -camtarget_Z]); camup([0 1 0]); axis off equal;
            set(gca,'CameraViewAngleMode', 'manual', 'CameraViewAngle', rad2deg(atan2(min(Xp_half_size),camera_calib.camera_array{view}.focal_length))*2);
            campos([0 0 -camtarget_Z]+[sind(0) 0 cosd(0)]*camtarget_Z); 
            camlight('headlight');

            msubplot(3,5,13,hs,vs,tb,bb,lb,rb);
            hold on;
            for view=1 %:size(ProjectionMatrices_pix,3)
                draw_coordinate_system(CameraPositions(:,:,view), 50, 3);
                surface(reshape(XYZ_Source{view}(1,:),size(X)),reshape(XYZ_Source{view}(2,:),size(X)),reshape(XYZ_Source{view}(3,:),size(X)), ...
                    'FaceColor','texturemap','EdgeColor','none', 'CData',repmat(mat2gray(Xp_imgs{folder_indx,view}(:,:,k)'),[1 1 3]), 'CDataMapping', 'direct', 'FaceLighting', 'none');
            end
            for i=1:length(target_objects)
                vert = ApplyTransform4x4(all_polygons{i}.vertices, diag([1 -1 1 1])*RegTools.convertTransRotTo4x4(all_estimated_6Nxnum_frame{folder_indx}(:,i,k)));
                patch('Faces', all_polygons{i}.faces, 'Vertices', vert, 'FaceColor', object_colors(target_objects(i),:), 'EdgeColor', 'none', 'FaceAlpha', 1.0); %object_alphas(indx_in_polygon)); %, 'EdgeColor', 'None', 'CData', u_disc, 'FaceColor', 'interp');
            end
            set(gca, 'clim', [0 10]);
            camera_target = all_estimated_6Nxnum_frame{folder_indx}(1:3,1,1)'; %all_moved3Dlandmark_center(1,:); % + [0 200 0];
            camproj('perspective'); camtarget(camera_target); camup([0 1 0]); axis off equal;
            set(gca,'CameraViewAngleMode', 'manual', 'CameraViewAngle', 20);
            material dull;
            view_angle = 50; camera_distance = 800;
            campos(camera_target+[sind(view_angle) 0.0 cosd(view_angle)]*camera_distance); 
            camlight('headlight');
            
            msubplot(3,5,15,hs,vs,tb,bb,lb,rb);
            hold on;
            base_bone_indx = 2;
            base_bone_transform = RegTools.convertTransRotTo4x4(all_estimated_6Nxnum_frame{folder_indx}(:,base_bone_indx,k));
            for i=1:length(target_objects)
                transform = RegTools.convertTransRotTo4x4(all_estimated_6Nxnum_frame{folder_indx}(:,i,k));
                vert = ApplyTransform4x4(all_polygons{i}.vertices, RegTools.matTranslation(-mean(all_polygons{base_bone_indx}.vertices))*inv(base_bone_transform)*transform);
                faces = all_polygons{i}.faces;
                if(i==1), faces( mean([vert(faces(:,1),3), vert(faces(:,2),3), vert(faces(:,3),3)],2) > 50,: ) = []; end
                patch('Faces', faces, 'Vertices', vert, 'FaceColor', object_colors(target_objects(i),:), 'EdgeColor', 'none', 'FaceAlpha', 1.0); %object_alphas(indx_in_polygon)); %, 'EdgeColor', 'None', 'CData', u_disc, 'FaceColor', 'interp');
            end
            set(gca, 'clim', [0 10]);
            camera_target = [0 0 0]; %all_moved3Dlandmark_center(1,:) + [200 0 0];
            camproj('perspective'); camtarget(camera_target); camup([0 0 1]); axis off equal;
            set(gca,'CameraViewAngleMode', 'manual', 'CameraViewAngle', 20);
            material dull;
            view_angle = 180; camera_distance = 300;
            campos(camera_target+[cosd(view_angle) sind(view_angle) 0]*camera_distance); 
            camlight('headlight');
end
%             colorbar('Position', [0.65 0.05 0.01 0.3],'FontSize',14);
if(0)

            msubplot(2,4,8,hs,vs,tb,bb,lb,rb);
            ylim = [-35 10];
%             ylim = [-50 50];
            plot(l_hip_trans_angle_stats');
            set(gca,'xlim',[1 num_frames],'ylim',ylim,'FontSize',14);
            ylabel({'relative translation, rotation of left femur', '[mm or deg]'}, 'FontSize', 12);
            line([k k], ylim, 'Color', 'k', 'LineWidth', 1);
            legend( sprintf('Tx (%.1f)',l_hip_trans_angle_stats_init(1,folder_indx)), sprintf('Ty (%.1f)',l_hip_trans_angle_stats_init(2,folder_indx)), ...
                    sprintf('Tz (%.1f)',l_hip_trans_angle_stats_init(3,folder_indx)), sprintf('Rx (%.1f)',l_hip_trans_angle_stats_init(4,folder_indx)), ...
                    sprintf('Ry (%.1f)',l_hip_trans_angle_stats_init(5,folder_indx)), sprintf('Rz (%.1f)',l_hip_trans_angle_stats_init(6,folder_indx)), ...
                    'Location', 'SouthEast');
end

%             btitle(sprintf('%s, %s, \\lambda_{L2}: %1.1e, Frame %d/%d', patient_ID, xray_image_files{1}, Lambda_L2_penalty, k, num_frames), 24, [0 0 0], 'tex', 0.97);
            btitle(sprintf('%s, %s & %s, Frame %d/%d', patient_ID, xray_image_files{:}, k, num_frames), 16, [0 0 0], 'tex', 0.97);
            drawnow;
            frame = getframe(gcf);
            writeVideo(v, frame);
            if(0) % save each frame as an image
                saveas(gcf, fullfile(output_rendering_folder_full, sprintf('frame_%d.png', k)));
            end
        end
        close(v);
    end
    
end