classdef CostParameters < hgsetget
    % Parameters used in ComputeCost function
    % This class should be 'light' since we pass it around many times in
    % the optimization process
    properties (Access = public)
        landmark_tbls = [];
%         num_landmarks_per_bone = [];
        ProjectionMatrices_pix = [];
        ctrl_scale = [];
        offset_scale = [];
        global_transform_6xMxN = [];
        num_ctrl_pnts = 1;
        num_transforms = 1;
        DimSize_2D = [512 512; 512 512];
        view_direction = [];
        landmark_tbls_2D_linear = [];
        landmark_tbls_3D_4v = [];
        enable_image = false;
        regTools = [];
        volumePlans = [];
        similarity_measure_plan_id = [];
        SimilarityMeasureType = 0;
        image_similarity_frame = [];
        volume_offset_4x4xN = [];
        rotation_center_local_4x4xN = [];
        surface_patch = [];
        surface_rendering_center_coordinate = [];
        surface_rendering_center_reference_object = [];
        previous_elapsed_time = 0;
        previous_iteration = 0;
        landmark_centroids_2D = [];
        lambda_rigidity_regularization = 1e-1;
        rigidity_constraint_object = [];
        lambda_smoothness_regularization = 0; %1e-1;
        joint_constraint_tbl = [];
        lambda_joint_regularization = 0; %1e-1;
        zoom_ROI_size = 50;
        maxParallelRendering = 20;
        GPU_IDs = [];
        CPU_par = false;
        PopSize = 10;
        fix_ctrl_offset = true;
        registration_frame = [];
        target_object = [];
        LCN_sigma = 0;
        landmark_cost_threshold = 0;
        cost_alpha = 0.5;
        cost_similarity_scaling = 1e1;
        cost_landmark_scaling = 1e-1;
        
        registration_results_folder = [];
        DS_2D = 1;
%         base_offset_param = [];
        
        initial_param = [];
        iteration_count = 0;
        previous_show_iteration = 0;
        progress_show_step = 1000;
        enable_output_debug_renderings = false;
        cost_log = [];
        current_layer = 0;
%         similarity_log = [];
        TimerID = [];
        writerObj = [];
        
        subject_ID = '';
        experiment_ID = '';
        
        % parameters for visualization of progress
        % default, for original image
        %  for plastic bone model:  fixed_clim = [0.6 1.2]; grad_clim = [-0.001 0.001]; moving_clim = [0 10];
        fixed_clim = [0 1.5];  
        grad_clim = [-0.005 0.005]; 
        moving_clim = [0 20];
        
% parameters for calibrationn
        pos_2D = [];
        pos_3D = [];
        focal_length_pix = 1000;
        origin_2D_pix = 100;
%         ProjectionMatrix_pix = [];
        CameraPosOffset = repmat(eye(4),[1 1 2]); % camera position offset (pre and post multiplication)
        calib_box_positions = [];
%         geomID = [];
%         
        img2D_ElementSpacing = [1 1];
%         img3D_ElementSpacing = [1 1 1];
%         img3D_planIDs = [];
%         
        bounding_box_3D = [0 0 0];
        param_scaling = [];
        RegistrationLogs = [];
        landmark_errors = [];
        
    end
end
