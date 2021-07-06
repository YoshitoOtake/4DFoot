addpath(getenv('RegToolsPath'));
addpath('../');
addpath('../util');
addpath('../debug_tools');
GPU_IDs = [0 1];

enable_SuperSloMo_interpolation = true;
% bone_mode = 'all_tarsal_bones';
% bone_mode = 'all_metatarsal_bones';
% bone_mode = 'tarsal_bones_distal';
bone_mode_array = { 'tarsal_bones_proximal',  'tarsal_bones_distal', 'metatarsal_bones', 'tibia', 'talus', 'calcaneus', 'navicular'};

dataset_ID_array = {'auto', 'manual', 'seg_auto_landmark_manual', 'seg_manual_landmark_auto', 'bone_model_manual', 'bone_model_RSA'};
timestamp = datestr(now,'yyyymmdd_HHMMSS');
figure('Position',[100 150 [2000 900]*10/10], 'PaperPositionMode', 'auto', 'Color', 'w');
if(~isempty(gcp('nocreate'))), delete(gcp('nocreate')); end     % parallel pool seems to conflict with RegTools, so we delete the pool every time

lambda_rigidity_regularization = 0; %1e-3;
lambda_joint_regularization = 0; %1e-2;

registration_results_root_dir_prefix = 'D:/Collaboration/Nara Medical University/Foot2D3D_registration_results/20210305_rigidity0_joint0_alpha05_ver1_Foot2D3D_registration_results';
for patient_indx = 1 %[1 2 3 4 5] %[-1 1 2 3 4 5] 
%     if(patient_indx<0), dataset_indx_array = 5; bone_mode_indx_try = 1; else, dataset_indx_array = [1 2 3 4]; bone_mode_indx_try = [1 2 3]; end
    if(patient_indx<0), dataset_indx_array = 5; bone_mode_indx_try = 1; else, dataset_indx_array = 2; bone_mode_indx_try = [1]; end
%     dataset_indx_array = [2]; bone_mode_indx_try = [4 5 6 7]; %[1 2 3];
    for dataset_indx = dataset_indx_array
        dataset_ID = dataset_ID_array{dataset_indx};
        setup_dataset;
        registration_results_root_dir = sprintf('%s_%s', registration_results_root_dir_prefix, dataset_ID);

        for bone_mode_indx = bone_mode_indx_try
            bone_mode = bone_mode_array{bone_mode_indx};
            fprintf('start footXrayVideo_NaraMedUniv, dataset_ID: %s, bone: %s\n', dataset_ID, bone_mode);
            bone_tbl = generate_bone_table( bone_mode );
            if(dataset_indx_array == 6), registration_mode_array = [1 3];
            else, if(bone_mode_indx ==1 || bone_mode_indx == 4), registration_mode_array = [1 2]; else, registration_mode_array = 2; end
            end
            for registration_mode = registration_mode_array %[1 2]  % 1: initial global registration, 2: fine local (frame-wise) registration, 3: fine local (frame-wise) registration without images for RSA
                [patient_ID, all_start_frames, all_end_frames, all_xray_image_files, calibration_patient_ID] = patient_specific_setup(patient_indx, enable_SuperSloMo_interpolation);

                for experiment_ID = 1 %:length(all_xray_image_files)
                    xray_image_files = all_xray_image_files{experiment_ID};
                    start_frame = all_start_frames(experiment_ID);
                    end_frame = all_end_frames(experiment_ID);
                    switch registration_mode
                        case 1
                            fprintf('run global registration\n');
                            global_registration_result_file = []; enable_image = false; fix_ctrl_offset = false; num_ctrl_pnts = 10;
                            registration_results_folder_mode_prefix = 'global';
                            bone_mode_prefix = '';
                        case {2, 3}
                            fprintf('run local registration\n');
                            global_registration_result_dir = getLatestTimestampedFile( fullfile(registration_results_root_dir, sprintf('results_global_%s_%s_%s_*',patient_ID,xray_image_files{1},xray_image_files{2})) );
                            global_registration_result_file = fullfile(global_registration_result_dir.folder, global_registration_result_dir.name, 'final_results.json');
                            enable_image = true; fix_ctrl_offset = true; 
                            if(registration_mode==2), registration_results_folder_mode_prefix = 'local'; else, registration_results_folder_mode_prefix = 'RSA_local'; end
                            bone_mode_prefix = ['_' bone_mode];
                    end
                    registration_results_folder = fullfile(registration_results_root_dir, sprintf('results_%s_%s_%s_%s%s_%s',registration_results_folder_mode_prefix,patient_ID,xray_image_files{1},xray_image_files{2},bone_mode_prefix,timestamp) );
                    if(~exist(registration_results_folder,'dir')), mkdir(registration_results_folder); end
                    if(~isempty(global_registration_result_file) && ~exist(global_registration_result_file, 'file')), fprintf('%s not found\n', global_registration_result_file); continue; end
                    fprintf('start registration on %s and %s, registration_mode: %d\n', xray_image_files{:}, registration_mode);

                    % load CT
                    if(enable_image)
                        mask_filter_width = 0; %4; %3;
                        [myu_img, HU_header, mask_img, HU_filename, mask_filename] = Load3DImage_foot(data_root_dir_image, patient_ID, CT_file_name);
                    else
                        HU_filename = []; mask_filename = [];
                        [~, HU_header] = mhdread(fullfile(data_root_dir_image, 'image', 'CT_left_leg_FC05', patient_ID, CT_file_name{1}, CT_file_name{2}) , true );  % load header only
                    end

                    % load Xray images
                    if(enable_image)
                        label_2D_erosion = -1; % 5; %15;
                        disable_log_correction = false;
                        border_mask = 5;
                        if(dataset_indx==5), med_filt_size = []; else, med_filt_size = [3 3]; end
                        [img_2D_original, similarity_measure_mask, lateral_filename, oblique_filename] = Load2DImages_foot(data_root_dir_image, patient_ID, Xray_suffix, view_direction, strcat(xray_image_files, x_ray_image_file_suffix), leading_view, label_2D_erosion, disable_log_correction, border_mask, med_filt_size);
                        LCN_sigma = 0; %0;
                        img_2D_original = ApplyLCN_image_wise(img_2D_original, LCN_sigma);
                    else
                        lateral_filename = [];
                        oblique_filename = [];
                    end

                    % load landmark files
                    [landmark_tbls, joint_center_tbl, landmark_2D_filename, landmark_3D_filename] = LoadLandmarkFiles_foot(data_root_dir_landmark, x_ray_landmark_dir, patient_ID, view_direction, xray_image_files, x_ray_landmark_filename, bone_tbl.Row, leading_view, CT_file_name, CT_landmark_filename);
                    if(dataset_indx==5),
                        % remove the second beads in tibia (seems to be moved between 2D and 3D scan)
                        for view=1:2, landmark_tbls{view}.pos3D{1}(2,:) = [];  landmark_tbls{view}.pos2D{1}(:,:,2) = []; end
                    end
                    num_views = size(landmark_tbls,1);
                    
                    % save filenames
                    S = struct('HU_filename', HU_filename, 'mask_filename', mask_filename, 'Xray_lateral_filename', lateral_filename, 'Xray_oblique_filename', oblique_filename, ...
                        'landmark_2D_lateral_filename', landmark_2D_filename{1}, 'landmark_2D_oblique_filename', landmark_2D_filename{2}, ...
                        'landmark_3D_lateral_filename', landmark_3D_filename{1}, 'landmark_3D_oblique_filename', landmark_3D_filename{2} ...
                        );
                    fid = fopen( fullfile(registration_results_folder, 'filenames.json'), 'w' );
                    fprintf(fid, '%s', prettyjson(jsonencode(S)));
                    fclose(fid);

                    % load calibration file
                    [ProjectionMatrices_pix_3D, DimSize_2D, camera_calib] = LoadCameraCalibration_foot(data_root_dir_image, calibration_patient_ID);

                    % preparation for optimization
                    opt_data = CostParameters();
                    opt_data.ProjectionMatrices_pix = cat(1, ProjectionMatrices_pix_3D(:,:,1), ProjectionMatrices_pix_3D(:,:,2));
                    opt_data.DimSize_2D = DimSize_2D;
                    opt_data.view_direction = view_direction;
                    opt_data.enable_image = enable_image;
                    opt_data.GPU_IDs = GPU_IDs;
                    opt_data.CPU_par = true;    % 1) parfor may crash Matlab in the use with RegTools, 2) multi CPU threads may not be always faster due to parallelization overhead
                    opt_data.fix_ctrl_offset = fix_ctrl_offset;
                    opt_data.registration_results_folder = registration_results_folder;
                    opt_data.progress_show_step = 50000; %500;
                    opt_data.SimilarityMeasureType = RegTools.SimilarityMeasureType_GC;
                    opt_data.fixed_clim = fixed_clim;
                    opt_data.grad_clim = grad_clim;
                    opt_data.moving_clim = moving_clim;

                    % initialize transformation parameters
                    num_total_frames = size(landmark_tbls{1}.pos2D{1},1);
                    if(isempty(global_registration_result_file))
                        % for global registration
    %                     opt_data.ctrl_scale = [1 1 1 1 1 1]*10;      % mm mm mm deg deg deg
    %                     opt_data.offset_scale = [1 1 1 1 1 1]*10;    % mm mm mm deg deg deg
                        opt_data.num_transforms = 1;
                        % centralize 3D landmark
                        pos3D_center_offset = [eye(3) -mean([cell2mat(landmark_tbls{1}.pos3D); cell2mat(landmark_tbls{2}.pos3D)])'; 0 0 0 1];
                        global_transform_4x4 = RegTools.convertTransRotTo4x4([0 0 -700 0 -90 0]) * pos3D_center_offset;
                        global_transform_6xMxN = repmat(RegTools.convert4x4ToTransRot(global_transform_4x4)', 1, num_total_frames, opt_data.num_transforms);
                    else
                        opt_data.num_transforms = height(bone_tbl); %height(landmark_tbls{1});  % one transform per bone
                        if( bone_mode_indx==2 )
                            ref_bone_indx = 3; 
                            global_transform_6xMxN = getGlobalTransformation_from_previous_registration( registration_results_root_dir, patient_ID, xray_image_files, bone_mode_array{bone_mode_indx-1}, ref_bone_indx, num_total_frames, opt_data.num_transforms);
                        elseif(bone_mode_indx==3)
                            ref_bone_indx = 1;
                            global_transform_6xMxN = getGlobalTransformation_from_previous_registration( registration_results_root_dir, patient_ID, xray_image_files, bone_mode_array{bone_mode_indx-1}, ref_bone_indx, num_total_frames, opt_data.num_transforms);
                        else
                            loaded = jsondecode(fileread(global_registration_result_file));
                            loaded.transformation_parameters = repmat(loaded.transformation_parameters, [1 1 opt_data.num_transforms 1]);
                            global_transform_6xMxN = permute(reshape(RegTools.convert4x4ToTransRot_multi(reshape(loaded.transformation_parameters,4,4,[])),6,opt_data.num_transforms,[]),[1 3 2]);
                        end
                        % rearrange landmark table
                        new_landmark_tbls = cell(num_views, opt_data.num_transforms);
                        for i=1:num_views
                            for j=1:opt_data.num_transforms
                                new_landmark_tbls{i,j} = landmark_tbls{i}(j,:);
                            end
                        end
                        landmark_tbls = new_landmark_tbls;
                    end

                    if(enable_image)
                        opt_data.LCN_sigma = LCN_sigma;
                        opt_data.regTools = RegTools(GPU_IDs, [], 'log_file.txt');
                    end

                    ROI_size = 150;
                    % select frames to register
                    switch registration_mode
                        case 1
            %         if(isempty(global_registration_result_file))
                            all_registration_frames = {1:num_total_frames};
                            DS = [1 1]; 
                            PopSize_array = 300;
                            maxParallelRendering_array = [];    % only for the case with image
                            MaxFunEvals = 1e6; %[1e4 3e4];
                            cost_alphas = 0.0;
                            search_range_array = 10;
                            opt_data.rigidity_constraint_object = [];
                            opt_data.lambda_rigidity_regularization = 0; %1e-1;
                            opt_data.lambda_smoothness_regularization = 0;
                        case 2
            %         else
                            num_registration_frames = 1; %5;
                            frame_stride = 1; %2;
                            num_ctrl_pnts = num_registration_frames;  % per frame control point
                            registration_start = max(1, start_frame-frame_stride);
                            registration_end = min(num_total_frames,end_frame);
                            all_registration_frames = cellfun(@(x) x:min(x+num_registration_frames-1,num_total_frames), num2cell(registration_start:frame_stride:registration_end)', 'UniformOutput', false);
                            % settings for multi-resolution pyramid
                            DS = [1 2; 1 1]; % downsample ratio in [3D 2D]
                            PopSize_array = [1000; 400];
%                             maxParallelRendering_array = [250 50]; % for num_registration_frames == 1
                            maxParallelRendering_array = [500 200]; %[50 20];
%                             if(patient_indx<0)
%                                 search_range_array = [1 0.5]*10;    % bone model (larger tibia-talus joint)
%                             else
%                                 search_range_array = [1 0.5]'*2;
                                search_range_array = [1 0.5]'*10;
%                             end
    %                         search_range_array = [search_range_array; search_range_array/2];
    %                         search_range_array = [ones(1,3)*0.01 ones(1,3)*70 ones(1,18)*0.01]; %[1 0.5]*5;
                            MaxFunEvals = [1e5 1e5]*2; %[1e4 3e4];
%                             cost_alphas = [0.1 1];
                            cost_alphas = [0.5 1];

                            % setting of rigidity constraints
%                             lambda_rigidity_regularization = 1e-3;
                            switch bone_mode
                                case {'tarsal_bones_proximal'}
                                    opt_data.rigidity_constraint_object = 2:4;
                                    opt_data.lambda_rigidity_regularization = lambda_rigidity_regularization;
                                case {'tarsal_bones_distal'}
                                    opt_data.rigidity_constraint_object = 1:4;
                                    opt_data.lambda_rigidity_regularization = lambda_rigidity_regularization;
%                                     DS = DS(2,:); PopSize_array = PopSize_array(2); maxParallelRendering_array = maxParallelRendering_array(2); 
%                                     search_range_array = search_range_array(2,:); MaxFunEvals = MaxFunEvals(2); 
                                    search_range_array = search_range_array/2;
                                    cost_alphas = [1 1]; %cost_alphas(2);
                                    ROI_size = 100;
                                case {'metatarsal_bones'}
                                    opt_data.rigidity_constraint_object = 1:4;
                                    opt_data.lambda_rigidity_regularization = lambda_rigidity_regularization;
%                                     DS = DS(2,:); PopSize_array = PopSize_array(2); maxParallelRendering_array = maxParallelRendering_array(2); 
%                                     search_range_array = search_range_array(2,:); MaxFunEvals = MaxFunEvals(2); 
                                    search_range_array = search_range_array/2;
                                    cost_alphas = [1 1]; % cost_alphas(2);
                                    ROI_size = 120;
                                case {'all_tarsal_bones'}
                                    opt_data.rigidity_constraint_object = 2:5;
                                    opt_data.lambda_rigidity_regularization = 1e-4;  %1e-0;
                                case {'all_metatarsal_bones'}
                                    opt_data.rigidity_constraint_object = 2:8;
                                    opt_data.lambda_rigidity_regularization = 1e-4;  %1e-0;
                            end
                            opt_data.lambda_smoothness_regularization = 0; %1e-2;

                        case 3  % RSA (metallic beads)
                            num_registration_frames = 1;
                            all_registration_frames = cellfun(@(x) x:min(x+num_registration_frames-1,num_total_frames), num2cell(1:ceil(num_registration_frames/2):(num_total_frames-ceil(num_registration_frames/2)))', 'UniformOutput', false);
                            % settings for multi-resolution pyramid
                            DS = [1 1]; % downsample ratio in [3D 2D]
                            PopSize_array = 500;
                            maxParallelRendering_array = 50;
                            search_range_array = 20;
                            MaxFunEvals = 3e5;
                            cost_alphas = 0.0;  % landmark only for RSA
                            num_ctrl_pnts = num_registration_frames;  % per frame control point
                            opt_data.lambda_rigidity_regularization = 0; % disable rigidity regularization
                            opt_data.lambda_smoothness_regularization = 0; % disable smoothness regularization
                    end

                    % preparation for 3D volumes
                %     opt_data.volume_offset_4x4xN = zeros(4,4,height(bone_tbl));
                    if(enable_image)
                        [all_masked_volumes, all_ElementSpacing_volume, opt_data.volume_offset_4x4xN] = prepareBoneVolumes(myu_img, mask_img, bone_tbl.label_indx, mask_filter_width, HU_header.ElementSpacing, opt_data.regTools, DS);
                    else
                        pos3D_center_offset = [eye(3) -mean([cell2mat(landmark_tbls{1}.pos3D); cell2mat(landmark_tbls{2}.pos3D)])'; 0 0 0 1];
                        opt_data.volume_offset_4x4xN = repmat(pos3D_center_offset,[1 1 height(bone_tbl)]);
                    end

                    % setting of joint center & joint constraints
                    opt_data.rotation_center_local_4x4xN = repmat(eye(4),[1 1 size(opt_data.volume_offset_4x4xN,3)]); %opt_data.volume_offset_4x4xN;  % use the volume center as a default joint center
                    if(registration_mode>=2)    % enable joint regularization only for local registration
                        available_joint_centers = intersect(bone_tbl.Row, joint_center_tbl.Row);
                        for i=1:length(available_joint_centers)    % replace metatarsal joints with the joint center landmark
                            bone_indx = bone_tbl{available_joint_centers{i},1};
                            opt_data.rotation_center_local_4x4xN(:,:,bone_indx) = opt_data.volume_offset_4x4xN(:,:,bone_indx) \ [eye(3) joint_center_tbl{i,1}'; 0 0 0 1];
                        end
                        opt_data.joint_constraint_tbl = GenerateJointConstraintTable_footXrayVideo_NaraMedUniv( bone_mode, opt_data.volume_offset_4x4xN, joint_center_tbl, bone_tbl );
    %                     opt_data.joint_constraint_tbl = GenerateJointConstraintTable_subtalar_joint( opt_data.volume_offset_4x4xN, joint_center_tbl, bone_tbl );
                        opt_data.lambda_joint_regularization = lambda_joint_regularization;
                    else
                        opt_data.lambda_joint_regularization = 0;
                    end

                    % save optimization settings
                    S = struct('DS', DS, 'PopSize_array', PopSize_array, 'maxParallelRendering_array', maxParallelRendering_array, 'search_range_array', search_range_array, ...
                        'MaxFunEvals', MaxFunEvals, 'cost_alphas', cost_alphas, 'lambda_smoothness_regularization', opt_data.lambda_smoothness_regularization, ...
                        'lambda_rigidity_regularization', opt_data.lambda_rigidity_regularization, 'rigidity_constraint_object', opt_data.rigidity_constraint_object, 'lambda_joint_regularization', opt_data.lambda_joint_regularization);
                    fid = fopen( fullfile(registration_results_folder, 'optimization_settings.json'), 'w' );
                    fprintf(fid, '%s', prettyjson(jsonencode(S)));
                    fclose(fid);

                    % start registration ("trial" refers to the frame number for local per frame registration)
                    for trial = 1:length(all_registration_frames)
                        opt_data.registration_frame = all_registration_frames{trial};
                        final_similarity_per_frame = [];
                        for target_object = 0 %1:height(bone_tbl)
                            if(isempty(global_registration_result_file))
                                title_suffix_str = '';
                                opt_data.target_object = 1;
                            else
                                if(target_object>0)
                                    opt_data.target_object = target_object; % target a specific bone
                                    title_suffix_str = sprintf('_%03d-%03d_%s', opt_data.registration_frame([1 end]), bone_tbl.Row{opt_data.target_object});
                                else
                                    opt_data.target_object = 1:height(bone_tbl); % target all bones;
                                    title_suffix_str = sprintf('_%03d-%03d_all_bones', opt_data.registration_frame([1 end]));
                                end
                            end
                            % initialization for one trial
                            opt_data.experiment_ID = sprintf('Patient:%s, Lateral:%s, Frame %d', patient_ID, xray_image_files{1}, trial);
                            opt_data.num_ctrl_pnts = num_ctrl_pnts; %length(opt_data.registration_frame);
                            if(enable_image)
                                img_2D = img_2D_original(:,:,:,opt_data.registration_frame);
                            end
                            opt_data.landmark_tbls = landmark_tbls;
                            for i=1:numel(landmark_tbls)
                                opt_data.landmark_tbls{i}.pos2D = cellfun(@(x) x(opt_data.registration_frame,:,:),landmark_tbls{i}.pos2D,'UniformOutput',false);
                            end
                            % calculate centroid of 2D landmarks at each frame
                            opt_data.landmark_tbls = opt_data.landmark_tbls(:,opt_data.target_object);
                            landmark_centroids_2D = CalculateLandmarkCentroids_per_frame(landmark_tbls, opt_data.target_object, global_transform_6xMxN, opt_data.volume_offset_4x4xN, ProjectionMatrices_pix_3D);

                            opt_data.global_transform_6xMxN = global_transform_6xMxN(:,opt_data.registration_frame,:);
                            opt_data.iteration_count = 0;
                            opt_data.previous_show_iteration = 0;
                            opt_data.cost_log = [];
                            opt_data.image_similarity_frame = 1:length(opt_data.registration_frame); 

                            opt_data.num_transforms = size(opt_data.landmark_tbls,2);
                            if(opt_data.fix_ctrl_offset)
                                opt_data.initial_param = zeros(opt_data.num_transforms*6*(opt_data.num_ctrl_pnts),1);
                            else
                                opt_data.initial_param = zeros(opt_data.num_transforms*6*(opt_data.num_ctrl_pnts+1),1);
                            end

                            visualization_on = true;
                            recording = false;
                            if(visualization_on)
                                if(recording)
                                    opt_data.writerObj = VideoWriter( fullfile(registration_results_folder, sprintf('optimization_process%s',title_suffix_str)) );
                                    opt_data.writerObj.FrameRate = 5;
                                    open(opt_data.writerObj);
                                else
                                    opt_data.writerObj = [];
                                end
                            end
                            opt_data.TimerID = tic;

                            for resolution_level=1:size(DS,1)
                                % run optimization
                                cmaes_opt = struct('MaxFunEvals', MaxFunEvals(resolution_level), 'TolFun', 1e-2, 'DispModulo', 0, 'LogModulo', 0, 'LogTime', 0, 'SaveVariables', 'off', ...
                                    'DispFinal', 'off', 'LogPlot', 'off', 'ReadSignals', 0, 'EvalInitialX', 0, 'EvalParallel', 'on', 'PopSize', 300);

                %                 opt_data.ctrl_scale = [1 1 1 5 5 5]*search_range_array(resolution_level);
                %                 opt_data.offset_scale = [1 1 1 5 5 5]*search_range_array(resolution_level);
                                if(numel(search_range_array)==size(DS,1))   % mm mm mm deg deg deg
                                    % one scalar value per resolution level
                                    opt_data.ctrl_scale = [1 1 1 1 1 1]*search_range_array(resolution_level);
                                    opt_data.offset_scale = [1 1 1 1 1 1]*search_range_array(resolution_level);
                                else
                                    % per each bone (6 x number of bones per resolution level)
                                    opt_data.ctrl_scale = search_range_array(resolution_level,:);
                                    opt_data.offset_scale = search_range_array(resolution_level,:);
                                end
                                if(enable_image)
                                    numGPUs = length(GPU_IDs);
                                    opt_data.DS_2D = DS(resolution_level,2);
                                    opt_data.maxParallelRendering = numGPUs * maxParallelRendering_array(resolution_level); %10;
                                    cmaes_opt.PopSize = floor(PopSize_array(resolution_level)/opt_data.maxParallelRendering)*opt_data.maxParallelRendering;    % PopSize needs to be a multiple of the number of GPUs
                                    num_image_sets = min(opt_data.maxParallelRendering, cmaes_opt.PopSize);

                                    fprintf('RegTools setup start: %.2f MB available on the GPU\n', opt_data.regTools.GPUmemCheck/1024/1024); timerID = tic();
                                    % prepare 2D images
                                    img_2D_DS = opt_data.regTools.Downsample2DProjections(DS(resolution_level,2), reshape(img_2D(:,:,:,opt_data.image_similarity_frame),size(img_2D,1),size(img_2D,2),[]));

                                    % prepare geometry
                                    geomID = opt_data.regTools.GenerateGeometry_3x4ProjectionMatrix( repmat(ProjectionMatrices_pix_3D, [1 1 num_image_sets*length(opt_data.image_similarity_frame)]), ...
                                        [1 1], [size(img_2D_DS,1) size(img_2D_DS,2)], [1 1]*DS(resolution_level,2) );

                                    % prepare for volume plan
                                    opt_data.volumePlans = zeros(height(bone_tbl),1,'int32');
                                    for j=1:height(bone_tbl)
                                        opt_data.volumePlans(j) = opt_data.regTools.CreateForwardProjectionPlan( all_masked_volumes{resolution_level,j}, all_ElementSpacing_volume(resolution_level,:) );
                                    end

                                    % prepare for similarity measure computation plan
                                    GI_Sigma = 1.0;
                                    mask = opt_data.regTools.Downsample2DProjections(DS(resolution_level,2), reshape(similarity_measure_mask(:,:,:,opt_data.registration_frame),size(img_2D,1),size(img_2D,2),[]));
                                    opt_data.similarity_measure_plan_id = opt_data.regTools.CreateSimilarityMeasureComputationPlan( img_2D_DS, GI_Sigma, mask, num_image_sets, [], 0, 0, 0, 0, [], opt_data.SimilarityMeasureType );
                                    fprintf('RegTools setup finished: %.2f MB available on the GPU\n', opt_data.regTools.GPUmemCheck/1024/1024);
                                    fprintf('initialization took %f sec\n', toc(timerID));

                                    % scale landmark centroids to fit current resolution level
                                    opt_data.landmark_centroids_2D = cellfun(@(x) x(opt_data.registration_frame,:)/DS(resolution_level,2),landmark_centroids_2D,'UniformOutput',false);
                                    opt_data.zoom_ROI_size = ROI_size/DS(resolution_level,2);
                                end
                                opt_data.PopSize = cmaes_opt.PopSize;
                                opt_data.cost_alpha = cost_alphas(resolution_level);

                                % show initial position
                                cost = landmark_intensity_registration_cost(opt_data.initial_param, opt_data, visualization_on, true);
                                if(0)
                                    % no optimization for debugging
                                    x_fitted = opt_data.initial_param;
                                else
                                    [x_fitted, cost_fitted, counteval] = ...
                                        cmaes_ex( 'landmark_intensity_registration_cost', ...
                                                opt_data.initial_param, ...  % objective variables initial point, determines N
                                                ones(length(opt_data.initial_param),1), ... % initial coordinate wise standard deviation(s)
                                                cmaes_opt, ...  % options struct, see defopts below
                                                opt_data, visualization_on);
                                end

                                % show optimized position and set the result to global_transform for next resolution level
                    %             enable_output_debug_renderings = (resolution_level==size(DS,1));
                                opt_data.enable_output_debug_renderings = false;
                                [cost, trans_param_4x4_image_coord] = landmark_intensity_registration_cost(x_fitted, opt_data, visualization_on, true);
                                opt_data.global_transform_6xMxN = permute(reshape(RegTools.convert4x4ToTransRot_multi(reshape(trans_param_4x4_image_coord,4,4,[])),6,[],length(opt_data.registration_frame)),[1 3 2]);
                                opt_data.initial_param = zeros(length(opt_data.initial_param),1);

                                if(enable_image)
                                    if(resolution_level == size(DS,1))
                                        [final_similarity_per_frame, all_floating_images] = calculate_similarity_per_frame(opt_data, similarity_measure_mask, size(img_2D), DS(resolution_level,2), img_2D_DS, opt_data.initial_param);
                                    end
                                    % cleanup regTools
                                    for j=1:height(bone_tbl), opt_data.regTools.DeleteForwardProjectionPlan( opt_data.volumePlans(j) ); end
                                    opt_data.regTools.DeleteProjectionParametersArray( geomID );
                                    opt_data.regTools.DeleteSimilarityMeasureComputationPlan( opt_data.similarity_measure_plan_id );
                                end

                                if(resolution_level == size(DS,1))
                                    print(gcf, '-r150', '-dpng', fullfile(registration_results_folder, sprintf('resolution_level_%d%s.png', resolution_level, title_suffix_str)));
%                                     if(enable_image)
%                                         mhdwrite(all_floating_images, fullfile(registration_results_folder, sprintf('floating_images_resolution_level_%d%s.mhd', resolution_level, title_suffix_str)), struct('ElementSpacing', [camera_calib.camera1.element_spacing' 1],'CompressedData',true));
%                                     end
                                end
                            end

                            if(~isempty(opt_data.writerObj))
                                close(opt_data.writerObj);
                            end

                            % output results to json file
                            S = struct('registration_frame', opt_data.registration_frame, 'num_ctrl_pnts', opt_data.num_ctrl_pnts, 'num_transforms', opt_data.num_transforms, 'num_frames', length(opt_data.registration_frame), ...
                                'ctrl_scale', opt_data.ctrl_scale, 'offset_scale', opt_data.offset_scale, ...
                                'volume_offset_4x4xN', opt_data.volume_offset_4x4xN, 'fix_ctrl_offset', opt_data.fix_ctrl_offset, ...
                                'global_transform_6xMxN', opt_data.global_transform_6xMxN, 'optimized_param', x_fitted, 'transformation_parameters', trans_param_4x4_image_coord, 'final_similarity_per_frame', final_similarity_per_frame);

                            fid = fopen( fullfile(registration_results_folder, sprintf('final_results%s.json', title_suffix_str)), 'w' );
                            fprintf(fid, '%s', prettyjson(jsonencode(S)));
                            fclose(fid);
                        end
                    end

                    if(enable_image)
                        if(opt_data.CPU_par) % parallel pool seems to conflict with RegTools, so we delete the pool every time
                            delete(gcp('nocreate'));
                        end
                        opt_data.regTools = [];
                        clear opt_data.regTools
                    end
                end
            end
        end
    end
end