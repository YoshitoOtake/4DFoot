dataset_ID = 'manual';
setup_dataset;
view_direction = {'lateral', 'oblique'};
leading_view = 2; % the "leading_view" starts earlier than the other view by 1 frame

label_2D_erosion = -1; % 5; %15;
disable_log_correction = false;
border_mask = 5;
all_img_2D = cell(5,1);
all_start_end = zeros(5,2);
for patient_indx = 1:5
    patient_ID = sprintf('ND%04d',patient_indx);
    fprintf('loading %s\n', patient_ID);
    switch patient_indx
        case 1, all_xray_image_files = { {'LL0027', 'LL0008'} }; calibration_patient_ID = patient_ID; start_frame = 1*2; end_frame = 19*2;
        case 2, all_xray_image_files = { {'LL0048', 'LL0020'} }; calibration_patient_ID = 'ND0002-5'; start_frame = 3*2; end_frame = 21*2;
        case 3, all_xray_image_files = { {'LL0022', 'LL0006'} }; calibration_patient_ID = 'ND0002-5'; start_frame = 4*2; end_frame = 19*2;
        case 4, all_xray_image_files = { {'LL0023', 'LL0007'} }; calibration_patient_ID = 'ND0002-5'; start_frame = 11*2; end_frame = 20*2;
        case 5, all_xray_image_files = { {'LL0016', 'LL0004'} }; calibration_patient_ID = 'ND0002-5'; start_frame = 12*2; end_frame = 24*2;
    end
    all_img_2D{patient_indx} = Load2DImages_foot(data_root_dir_image, patient_ID, view_direction, all_xray_image_files{1}, leading_view, label_2D_erosion, disable_log_correction, border_mask);
    all_start_end(patient_indx,:) = [start_frame end_frame];
end

%%
figure('Position',[100 150 [1600 900]*10/10], 'PaperPositionMode', 'auto', 'Color', 'w'); colormap(gray(256));
fixed_clim = [0 1.5];
for patient_indx = 1:5
    for view = 1:2
        msubplot(5,2,(patient_indx-1)*2+view);
        frame_indx = floor(linspace(all_start_end(patient_indx,1),all_start_end(patient_indx,2),5));
        im(permute(squeeze(all_img_2D{patient_indx}(:,:,view,frame_indx)),[2 1 3]), 'nrows', 1, 'labelslice', frame_indx);
        set(gca,'clim',fixed_clim);
        axis off;
        title(sprintf('ND%04d, %s',patient_indx, view_direction{view}))
    end
end
