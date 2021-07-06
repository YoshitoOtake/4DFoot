function [myu_img, HU_header, mask_img, HU_filename, mask_filename] =  Load3DImage_foot(root_dir, patient_ID, CT_file_name)
    fprintf('loading 3D images of %s...', patient_ID); timerID = tic;
    % load CT (CT header is important for offset of 3D landmark files)
    HU_filename = fullfile(root_dir, 'image/CT_left_leg_FC05', patient_ID, CT_file_name{1}, CT_file_name{2});
    [HU_img, HU_header] = mhdread( HU_filename );
    myu_water = 0.2683;     % see RegTools.HU2Myu() function for detail
%     myu_img = HU2Myu_nonlinear(HU_img+00, myu_water, 0); %RegTools.HU2Myu(HU_img, myu_water);
    myu_img = HU2Myu(HU_img-00, myu_water); %RegTools.HU2Myu(HU_img, myu_water);
%     myu_img(myu_img<0) = 0;
    mask_filename = fullfile(root_dir, 'image/CT_left_leg_FC05', patient_ID, CT_file_name{1}, 'label.mhd');
    [mask_img, ~] = mhdread( mask_filename );
    mask_img = removeLabelNoise(mask_img);
%     mask_img(1:floor(size(mask_img,1)/2),:,:) = 0;  % remove right foot
    if(0)   % for debugging
        figure; im(flip(sum(permute(myu_img(257:end,:,:).*(mask_img(257:end,:,:)>0),[2,3,1]),3),2)); colorbar; colormap(gray(256)); set(gca,'clim',[0,10],'DataAspectRatio',HU_header.ElementSpacing([2 3 1]));
        figure; im(flip(sum(permute(myu_img(257:end,:,:),[2,3,1]),3),2)); colorbar; colormap(gray(256)); set(gca,'clim',[0,10],'DataAspectRatio',HU_header.ElementSpacing([2 3 1]));
        
        %%
        figure; msubplot(1,2,1);
        im_contour_overlay(permute(HU_img(:,:,200:10:300),[2 1 3]), permute(mask_img(:,:,200:10:300),[2 1 3]));
        msubplot(1,2,2);
        im_contour_overlay(permute(HU_img(:,:,200:10:300),[2 1 3]), permute(mask_img(:,:,200:10:300),[2 1 3]), [], 'fill_region', 1, 'opacity', 0.8);
    end
%     mask_img(:,:,450:end) = 0;  % remove higher part of tibia-fibula
    fprintf(' done in %f sec\n', toc(timerID));