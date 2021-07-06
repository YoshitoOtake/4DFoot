function xray_image_files = getAllXrayImageFiles(data_root_dir, patient_ID)

lateral_xray = dir( fullfile(data_root_dir, 'image', 'Xray_interpolated', patient_ID, 'lateral', '*.mhd') );
oblique_xray = dir( fullfile(data_root_dir, 'image', 'Xray_interpolated', patient_ID, 'oblique', '*.mhd') );
lateral_IDs = cellfun(@(x) x{2},regexp({lateral_xray.name},'_','split'),'UniformOutput',false);
oblique_IDs = cellfun(@(x) x{2},regexp({oblique_xray.name},'_','split'),'UniformOutput',false);
xray_image_files = cellfun(@(x,y) {x,y}, lateral_IDs, oblique_IDs, 'UniformOutput', false);
