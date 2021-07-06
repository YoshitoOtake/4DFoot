addpath(getenv('RegToolsPath'));
addpath('../');
addpath('../util');

for dataset_indx =1:3
    switch dataset_indx
        case 1
            % for ND0001
            input_root_dir = '\\scallop\User\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data';
            calibration_setting_filename = fullfile(input_root_dir, 'phantom_center', 'ND0001', 'calibration_setting.json');
            output_folders = {'\\scallop\User\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data\calibration_results', '\\scallop\User\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data_predicted\calibration_results'};
        case 2
            % for ND0002-5
            input_root_dir = '\\scallop\User\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data';
            calibration_setting_filename = fullfile(input_root_dir, 'phantom_center', 'ND0002-5', 'calibration_setting.json');
            output_folders = {'\\scallop\User\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data\calibration_results', '\\scallop\User\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data_predicted\calibration_results'};
        case 3
            % for bone model experiment
            input_root_dir = '\\scallop\User\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data_model';
            calibration_setting_filename = fullfile(input_root_dir, 'phantom_center', 'MODEL1', 'calibration_setting.json');
            output_folders = {'\\scallop\User\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data_model\calibration_results'};
    end

    setting = jsondecode( fileread(calibration_setting_filename) );
    input_root_dir = setting.root_dir;
    patient_ID = setting.patient_ID;
    phantom_beads_position = setting.phantom_beads_position;
    lateral_landmark_files = setting.lateral_landmark_files;
    oblique_landmark_files = setting.oblique_landmark_files;

    % image property obtained from DICOM header
    [img_lateral, hdr_lateral] = mhdread( fullfile(input_root_dir, 'phantom_image', patient_ID, sprintf('BP_%s.mhd', lateral_landmark_files{1})) );
    [img_oblique, hdr_oblique] = mhdread( fullfile(input_root_dir, 'phantom_image', patient_ID, sprintf('BP_%s.mhd', oblique_landmark_files{1})) );
    initial_focal_length_mm = [hdr_lateral.SDD hdr_oblique.SDD];
    DimSize_2D = hdr_lateral.DimSize(1:2);
    ElementSpacing = hdr_lateral.ElementSpacing(1:2); 
    % FieldOfViewDimensions = [305 305];
    % ElementSpacing = FieldOfViewDimensions./DimSize_2D;

    output_video_filename = 'two_view_calibration_run.mp4';

    % load landmark positions (detected by circle fitting of the metallic
    % sphere at the corners of the calibration box)
    all_pos_2D = cell(length(lateral_landmark_files),2);
    for i=1:length(lateral_landmark_files)
        all_pos_2D{i,1} = load_phantom_points_file( fullfile(input_root_dir, 'phantom_center', patient_ID, [lateral_landmark_files{i} '.csv']) );
        all_pos_2D{i,2} = load_phantom_points_file( fullfile(input_root_dir, 'phantom_center', patient_ID, [oblique_landmark_files{i} '.csv']) );
    end

    % preparation of parameters for optimization
    opt_data = CostParameters();
    opt_data.pos_3D =phantom_beads_position;

    opt_data.img2D_ElementSpacing = ElementSpacing;
    focal_length_pix = initial_focal_length_mm ./ mean(opt_data.img2D_ElementSpacing);   % we assume a fixed (known) focal length and same for X- and Y- axes
    opt_data.origin_2D_pix = DimSize_2D./2;
    init_SOD = mean([hdr_lateral.SOD; hdr_oblique.SOD]); %600;

    opt_data.progress_show_step = 30000;

    visualization_on = true;
    recording = true;
    if(visualization_on)
        figure('Position',[100 100 [1600 900]*8/10], 'PaperPositionMode', 'auto', 'Color', 'w');
        if(recording)
            opt_data.writerObj = VideoWriter(output_video_filename);
            opt_data.writerObj.FrameRate = 5;
            open(opt_data.writerObj);
        else
            opt_data.writerObj = [];
        end
    end
    opt_data.TimerID = tic;

    options = optimset();
    options.Display = 'iter';
    options.MaxFunEvals = 1e7;

    cmaes_opt = struct('MaxFunEvals', options.MaxFunEvals, 'TolFun', 1e-4, 'DispModulo', 0, 'LogModulo', 0, 'LogTime', 0, 'SaveVariables', 'off', ...
        'DispFinal', 'off', 'LogPlot', 'off', 'ReadSignals', 0, 'EvalInitialX', 0, 'EvalParallel', 'on', 'PopSize', 500);
    run_initialization = 1;
    % initialization
    if(run_initialization)
    %     all_initial_estimate = zeros([7 size(all_pos_2D)]);
        all_initial_estimate = zeros([6 size(all_pos_2D)]);
        for i=1:size(all_pos_2D,1)
            for j=1:size(all_pos_2D,2)
                opt_data.focal_length_pix(1) = focal_length_pix(j);
                opt_data.initial_param = [0 0 -init_SOD 0 0 -90];
                opt_data.param_scaling = ones(1,6);
                opt_data.pos_2D = all_pos_2D(i,j);
                opt_data.iteration_count = 0; opt_data.previous_show_iteration = 0; opt_data.cost_log = []; opt_data.RegistrationLogs=[];

                optimization_mode = 1;  % 1: 6-DOF, 2: 7-DOF (+focal length), 3: 
                % show initial position
                cost = two_view_calibration_cost(zeros(length(opt_data.initial_param),1), opt_data, visualization_on, optimization_mode);

                [x_fitted, cost_fitted, counteval] = ...
                    cmaes_ex( 'two_view_calibration_cost', ...
                            zeros(length(opt_data.initial_param),1), ...  % objective variables initial point, determines N
                            ones(length(opt_data.initial_param),1), ... % initial coordinate wise standard deviation(s)
                            cmaes_opt, ...  % options struct, see defopts below
                            opt_data, visualization_on, optimization_mode);

    %             opt_data.iteration_count = 0; opt_data.previous_show_iteration = 0; opt_data.cost_log = []; opt_data.RegistrationLogs=[];
    %             estimated = opt_data.initial_param(:) + x_fitted(:) .* opt_data.param_scaling(:);
    %             opt_data.initial_param = [focal_length_pix(j) estimated'];
    %             opt_data.param_scaling = ones(1,7);
    %             optimization_mode = 2;
    % 
    %             % show initial position
    %             cost = two_view_calibration_cost(zeros(length(opt_data.initial_param),1), opt_data, visualization_on, optimization_mode);
    % 
    %             [x_fitted, cost_fitted, counteval] = ...
    %                 cmaes_ex( 'two_view_calibration_cost', ...
    %                         zeros(length(opt_data.initial_param),1), ...  % objective variables initial point, determines N
    %                         ones(length(opt_data.initial_param),1), ... % initial coordinate wise standard deviation(s)
    %                         cmaes_opt, ...  % options struct, see defopts below
    %                         opt_data, visualization_on, optimization_mode);

                % show optimized position
                cost = two_view_calibration_cost(x_fitted, opt_data, visualization_on, optimization_mode, true);

                all_initial_estimate(:,i,j) = opt_data.initial_param(:) + x_fitted(:) .* opt_data.param_scaling(:);
    %             all_initial_estimate(4:6,i,j) = mod(all_initial_estimate(4:6,i,j), 360);    % clean angle parameter
            end
        end
        save('temp_initialization.mat', 'all_initial_estimate');
    else
        load('temp_initialization.mat');
    end

    % fine-tuning of relative camera position
    % focal_length_pix = squeeze(mean(all_initial_estimate(1,:,:),2))';
    % all_initial_estimate = all_initial_estimate(2:end,:,:);
    num_images = size(all_initial_estimate,2);
    all_relative_trans = zeros(6,num_images);
    for i=1:num_images
        all_relative_trans(:,i) = RegTools.convert4x4ToTransRot( ...
                RegTools.convertTransRotTo4x4(all_initial_estimate(:,i,1))*inv(RegTools.convertTransRotTo4x4(all_initial_estimate(:,i,2))) );
    end
    relative_pos_init = mean(all_relative_trans,2);    % simple mean of Euler angle may not work in some cases... (Rodriguez angle may be better)
    opt_data.CameraPosOffset(:,:,1) = RegTools.convertTransRotTo4x4( all_initial_estimate(:,1,1) );
    opt_data.CameraPosOffset(:,:,2) = opt_data.CameraPosOffset(:,:,1) \ RegTools.convertTransRotTo4x4( relative_pos_init );
    opt_data.focal_length_pix = focal_length_pix;
    opt_data.calib_box_positions = all_initial_estimate(:,:,1);
    opt_data.initial_param = [focal_length_pix DimSize_2D./2 DimSize_2D./2 zeros(1,6) reshape(all_initial_estimate(:,:,1),1,[])];
    opt_data.param_scaling = [ones(1,2)*10 ones(1,4)*10 ones(1,6) ones(1,6*num_images)*1];
    opt_data.pos_2D = all_pos_2D;
    opt_data.iteration_count = 0; opt_data.previous_show_iteration = 0; opt_data.cost_log = []; opt_data.RegistrationLogs=[];

    % show initial position
    % optimization_mode = 3; PM_offset = 0;
    optimization_mode = 4; PM_offset = 6;
    cost = two_view_calibration_cost(zeros(length(opt_data.initial_param),1), opt_data, visualization_on, optimization_mode);

    [x_fitted, cost_fitted, counteval] = ...
        cmaes_ex( 'two_view_calibration_cost', ...
                zeros(length(opt_data.initial_param),1), ...  % objective variables initial point, determines N
                ones(length(opt_data.initial_param),1), ... % initial coordinate wise standard deviation(s)
                cmaes_opt, ...  % options struct, see defopts below
                opt_data, visualization_on, optimization_mode);
    % disp(x_fitted');
    % show optimized position
    cost = two_view_calibration_cost(x_fitted, opt_data, visualization_on, optimization_mode, true);

    estimated_param = opt_data.initial_param(:) + x_fitted(:) .* opt_data.param_scaling(:);
    estimated_relative_pos_4x4 = opt_data.CameraPosOffset(:,:,1) * RegTools.convertTransRotTo4x4(estimated_param((PM_offset+1):end)) * opt_data.CameraPosOffset(:,:,2);
    if(optimization_mode == 4)
        focal_length_pix = estimated_param(1:2);
        origin_2D_pix_camera1 = estimated_param(3:4);
        origin_2D_pix_camera2 = estimated_param(5:6);
        opt_data.calib_box_positions = reshape(estimated_param((PM_offset+6+1):end),6,[]);
    else
        focal_length_pix = [opt_data.focal_length_pix(1), opt_data.focal_length_pix(2)];
        origin_2D_pix_camera1 = opt_data.origin_2D_pix;    
        origin_2D_pix_camera2 = opt_data.origin_2D_pix;
    end
    initial_focal_length_mm = focal_length_pix*mean(opt_data.img2D_ElementSpacing);
    fprintf('Estimated relative pos: %f, %f, %f, %f, %f, %f, focal_length_mm: %f, %f, origin: (%f, %f) , (%f, %f)\n', RegTools.convert4x4ToTransRot( estimated_relative_pos_4x4 ), initial_focal_length_mm, origin_2D_pix_camera1, origin_2D_pix_camera2);

    if(~isempty(opt_data.writerObj))
        close(opt_data.writerObj);
    end

    output_struct = struct( ...
        'focal_length',{initial_focal_length_mm(1), initial_focal_length_mm(2)}, ...
        'image_center',{origin_2D_pix_camera1, origin_2D_pix_camera2}, ...
        'final_landmark_error',{opt_data.landmark_errors(:,1), opt_data.landmark_errors(:,2)}, ...
        'element_spacing', {ElementSpacing, ElementSpacing}, ...
        'dim_size', {DimSize_2D, DimSize_2D} ...
        );
    json_text = jsonencode(struct('camera1', output_struct(1), 'camera2', output_struct(2), 'relative_position', estimated_relative_pos_4x4', 'calib_box_position',opt_data.calib_box_positions));

    for i=1:length(output_folders)
        output_json_folder = output_folders{i}; % fullfile(root_dir, 'calibration_results');
        if(~exist(output_json_folder,'dir')), mkdir(output_json_folder); end
        output_json_filename = fullfile(output_json_folder, sprintf('calibration_result_%s.json',patient_ID));

        fid = fopen(output_json_filename, 'w');
        fprintf(fid, prettyjson(json_text));
        fclose(fid);
        fprintf('saved calibration results to %s\n', output_json_filename);
    end
end

function pos_2D = load_phantom_points_file(filename)
    fid = fopen(filename);
    pos_2D = cell2mat(textscan(fid, '%f %f', 'Delimiter', ','));
    fclose(fid);
end
