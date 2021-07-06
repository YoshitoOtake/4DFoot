function [ProjectionMatrices_pix_3D, DimSize_2D, camera_calib] = LoadCameraCalibration_foot(root_dir, patient_ID)

camera_calib = jsondecode(fileread(fullfile(root_dir, 'calibration_results', sprintf('calibration_result_%s.json', patient_ID))));
adjustment_ratio1 = 1.0; % adjustment for either focal length or element spacing
adjustment_ratio2 = 1.0; % adjustment for either focal length or element spacing
camera_calib.camera1.focal_length = camera_calib.camera1.focal_length * adjustment_ratio1;
camera_calib.camera2.focal_length = camera_calib.camera2.focal_length * adjustment_ratio2;
f1 = (camera_calib.camera1.focal_length) / mean(camera_calib.camera1.element_spacing);
f2 = (camera_calib.camera2.focal_length) / mean(camera_calib.camera2.element_spacing);
ProjectionMatrix1_pix = [[-f1, 0, camera_calib.camera1.image_center(1); 0, -f1, camera_calib.camera1.image_center(2); 0, 0, 1] [0 0 0]'];
ProjectionMatrix2_pix = [[-f2, 0, camera_calib.camera2.image_center(1); 0, -f2, camera_calib.camera2.image_center(2); 0, 0, 1] [0 0 0]'] / (camera_calib.relative_position');
ProjectionMatrices_pix_3D = cat(3, ProjectionMatrix1_pix, ProjectionMatrix2_pix);
DimSize_2D = [camera_calib.camera1.dim_size(:)'; camera_calib.camera2.dim_size(:)'];
