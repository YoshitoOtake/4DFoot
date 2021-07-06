for dataset_indx = 1:3
    switch dataset_indx
        case 1
            root_dir = '\\scallop\User\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data\';
            patient_ID = 'ND0001';
            lateral_landmark_files = {'LL0020', 'LL0021', 'LL0022'};
            oblique_landmark_files = {'LL0001', 'LL0002', 'LL0003'};
            phantom_beads_position = [-1 -1 1; 1 -1 1; 1 1 1; -1 1 1; -1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1] * 55;
        case 2
            root_dir = '\\scallop\User\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data\';
            patient_ID = 'ND0002-5';
            lateral_landmark_files = {'LL0029', 'LL0030', 'LL0031', 'LL0032', 'LL0033', 'LL0034', 'LL0035', 'LL0036', 'LL0037', 'LL0038', 'LL0039', 'LL0040'};
            oblique_landmark_files = {'LL0001', 'LL0002', 'LL0003', 'LL0004', 'LL0005', 'LL0006', 'LL0007', 'LL0008', 'LL0009', 'LL0010', 'LL0011', 'LL0012'};
            phantom_beads_position = [-1 -1 1; 1 -1 1; 1 1 1; -1 1 1; -1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1] * 55;
        case 3
            root_dir = '\\scallop\User\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data_model\';
            patient_ID = 'MODEL1';
            lateral_landmark_files = {'LL0022', 'LL0023', 'LL0024', 'LL0025', 'LL0026', 'LL0027', 'LL0028', 'LL0029', 'LL0030', 'LL0031', 'LL0032', 'LL0033'};
            oblique_landmark_files = {'LL0002', 'LL0003', 'LL0004', 'LL0005', 'LL0006', 'LL0007', 'LL0008', 'LL0009', 'LL0010', 'LL0011', 'LL0012', 'LL0013'};
            phantom_beads_position = [-1 -1 1; 1 -1 1; 1 1 1; -1 1 1; -1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1] * 55;
            phantom_beads_position = phantom_beads_position([1 4 3 2 5 8 7 6],:); % mirror for bone model experiment
    end

    S = struct('root_dir', root_dir, 'patient_ID', patient_ID, 'phantom_beads_position', phantom_beads_position, 'lateral_landmark_files', {lateral_landmark_files}, 'oblique_landmark_files', {oblique_landmark_files});
    output_filename = fullfile( root_dir, 'phantom_center', patient_ID, 'calibration_setting.json');

    fid = fopen( output_filename, 'w' );
    fprintf(fid, '%s', prettyjson(jsonencode(S)));
    fclose(fid);
    fprintf('saved calibration setting in %s\n', output_filename);
end