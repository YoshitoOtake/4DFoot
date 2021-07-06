function [img_2D_original, similarity_measure_mask, lateral_filename, oblique_filename] = Load2DImages_foot(root_dir, patient_ID, Xray_suffix, view_direction, xray_image_files, leading_view, label_2D_erosion, disable_log_correction, border_mask, med_filt_size)
    if(~exist('disable_log_correction','var')), disable_log_correction = false; end
    lateral_filename = fullfile(root_dir, 'image', sprintf('Xray%s',Xray_suffix), patient_ID, view_direction{1}, ['BP_' xray_image_files{1} Xray_suffix '.mhd']);
    oblique_filename = fullfile(root_dir, 'image', sprintf('Xray%s',Xray_suffix), patient_ID, view_direction{2}, ['BP_' xray_image_files{2} Xray_suffix '.mhd']);
    [img_2D_lateral, ~] = mhdread( lateral_filename );
    [img_2D_oblique, ~] = mhdread( oblique_filename );
    if(label_2D_erosion>=0)
        [label_2D_lateral, ~] = mhdread( fullfile(root_dir, 'image/Xray_skin', patient_ID, view_direction{1}, ['BP_' xray_image_files{1}], ['skin' Xray_suffix '.mhd']) );
        [label_2D_oblique, ~] = mhdread( fullfile(root_dir, 'image/Xray_skin', patient_ID, view_direction{2}, ['BP_' xray_image_files{2}], ['skin' Xray_suffix '.mhd']) );
    end
    if(leading_view==0)
        % no modification
    elseif(leading_view==1)
        img_2D_lateral = img_2D_lateral(:,:,2:1:end);
        img_2D_oblique = img_2D_oblique(:,:,1:1:end-1);
        if(label_2D_erosion>=0)
            label_2D_lateral = label_2D_lateral(:,:,2:1:end);
            label_2D_oblique = label_2D_oblique(:,:,1:1:end-1);
        end
    else
        img_2D_lateral = img_2D_lateral(:,:,1:1:end-1);
        img_2D_oblique = img_2D_oblique(:,:,2:1:end);
        if(label_2D_erosion>=0)
            label_2D_lateral = label_2D_lateral(:,:,1:1:end-1);
            label_2D_oblique = label_2D_oblique(:,:,2:1:end);
        end
    end
    num_frames = min(size(img_2D_lateral,3), size(img_2D_oblique,3));
    img_2D_lateral = img_2D_lateral(:,:,1:num_frames);
    img_2D_oblique = img_2D_oblique(:,:,1:num_frames);
    if(label_2D_erosion>=0)
        label_2D_lateral = label_2D_lateral(:,:,1:num_frames);
        label_2D_oblique = label_2D_oblique(:,:,1:num_frames);
    end
    if(size(img_2D_lateral,3) ~= size(img_2D_oblique,3))
        fprintf('video frame does not match\n');
        return;
    end
    [nu, nv, num_frame] = size(img_2D_lateral);
    if(disable_log_correction)
        img_2D_original = cat(3, reshape(single(img_2D_lateral), [nu nv 1 num_frame]), reshape(single(img_2D_oblique), [nu nv 1 num_frame]));
    else
        img_2D_original = cat(3, reshape(LogCorrection(single(img_2D_lateral),255*ones(num_frame,1)), [nu nv 1 num_frame]), reshape(LogCorrection(single(img_2D_oblique),255*ones(num_frame,1)), [nu nv 1 num_frame]));
    end
    
    if(~isempty(med_filt_size))
        for i=1:num_frame
            img_2D_original(:,:,1,i) = medfilt2(img_2D_original(:,:,1,i), med_filt_size);
            img_2D_original(:,:,2,i) = medfilt2(img_2D_original(:,:,2,i), med_filt_size);
        end
    end
    
    if(label_2D_erosion>0)
    %     gauss_sigma = 2;
        for i=1:num_frame
            label_2D_lateral(:,:,i) = imerode(label_2D_lateral(:,:,i),strel('disk',label_2D_erosion,4));
            label_2D_oblique(:,:,i) = imerode(label_2D_oblique(:,:,i),strel('disk',label_2D_erosion,4));
    %         img_2D_original(:,:,1,i) = imgaussfilt(img_2D_original(:,:,1,i),gauss_sigma);
    %         img_2D_original(:,:,2,i) = imgaussfilt(img_2D_original(:,:,2,i),gauss_sigma);
        end
    end
    if(label_2D_erosion>=0)
        similarity_measure_mask = permute(cat(4, label_2D_lateral, label_2D_oblique),[1 2 4 3]);
    else
        similarity_measure_mask = ones(size(img_2D_original)); 
        if(border_mask>0)
            similarity_measure_mask(1:border_mask,:,:) = 0;
            similarity_measure_mask(end-border_mask:end,:,:) = 0;
            similarity_measure_mask(:,1:border_mask,:) = 0;
            similarity_measure_mask(:,end-border_mask:end,:) = 0;
        end
    end
%     similarity_measure_mask(370:end,:,1,:) = 0;   % mask of the floor in lateral view

