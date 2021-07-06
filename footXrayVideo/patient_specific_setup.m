function [patient_ID, start_frame, end_frame, all_xray_image_files, calibration_patient_ID] = patient_specific_setup(patient_indx, enable_SuperSloMo_interpolation)

if(~exist('enable_SuperSloMo_interpolation','var')), enable_SuperSloMo_interpolation = false; end

if(patient_indx>0)
    patient_ID = sprintf('ND%04d',patient_indx);
    if(enable_SuperSloMo_interpolation), start_end_frame_ratio = 2; else, start_end_frame_ratio = 1; end
    
    % for MICCAI2021 submission videos
%     switch patient_indx
%         case 1, all_xray_image_files = { {'LL0027', 'LL0008'} }; calibration_patient_ID = patient_ID; start_frame = 1*start_end_frame_ratio; end_frame = 19*start_end_frame_ratio;
%         case 2, all_xray_image_files = { {'LL0048', 'LL0020'} }; calibration_patient_ID = 'ND0002-5'; start_frame = 3*start_end_frame_ratio; end_frame = 21*start_end_frame_ratio;
%         case 3, all_xray_image_files = { {'LL0022', 'LL0006'} }; calibration_patient_ID = 'ND0002-5'; start_frame = 4*start_end_frame_ratio; end_frame = 19*start_end_frame_ratio;
%         case 4, all_xray_image_files = { {'LL0023', 'LL0007'} }; calibration_patient_ID = 'ND0002-5'; start_frame = 11*start_end_frame_ratio; end_frame = 20*start_end_frame_ratio;
%         case 5, all_xray_image_files = { {'LL0016', 'LL0004'} }; calibration_patient_ID = 'ND0002-5'; start_frame = 12*start_end_frame_ratio; end_frame = 24*start_end_frame_ratio;
%     end
    
    switch patient_indx
        case 1, all_xray_image_files = { {'LL0029', 'LL0010'}, {'LL0037', 'LL0019'} }; calibration_patient_ID = patient_ID; start_frame = [3 17]*start_end_frame_ratio; end_frame = [21 34]*start_end_frame_ratio;
        case 2, all_xray_image_files = { {'LL0048', 'LL0020'}, {'LL0056', 'LL0028'} }; calibration_patient_ID = 'ND0002-5'; start_frame = [3 7]*start_end_frame_ratio; end_frame = [20 19]*start_end_frame_ratio;
        case 3, all_xray_image_files = { {'LL0023', 'LL0007'}, {'LL0030', 'LL0014'} }; calibration_patient_ID = 'ND0002-5'; start_frame = [4 6]*start_end_frame_ratio; end_frame = [17 18]*start_end_frame_ratio;
        case 4, all_xray_image_files = { {'LL0024', 'LL0008'}, {'LL0032', 'LL0016'} }; calibration_patient_ID = 'ND0002-5'; start_frame = [11 10]*start_end_frame_ratio; end_frame = [23 19]*start_end_frame_ratio;
        case 5, all_xray_image_files = { {'LL0017', 'LL0005'}, {'LL0024', 'LL0012'} }; calibration_patient_ID = 'ND0002-5'; start_frame = [4 5]*start_end_frame_ratio; end_frame = [16 18]*start_end_frame_ratio;
    end

else
    patient_ID = 'MODEL1'; all_xray_image_files = { {'LL0038', 'LL0018'} };  calibration_patient_ID = patient_ID; start_frame = 1; end_frame = 87;
end
