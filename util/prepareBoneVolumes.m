function [all_masked_volumes, all_ElementSpacing_volume, volume_offset_4x4xN] = ...
    prepareBoneVolumes(myu_img, mask_img, bone_label_indx, mask_filter_width, ElementSpacing, regTools, DS)
    num_bones = length(bone_label_indx);
    volume_offset_4x4xN = zeros(4,4,num_bones);
    all_masked_volumes = cell(size(DS,1),num_bones);
    all_ElementSpacing_volume = zeros(size(DS,1),3);
    fprintf('preparing bone volumes'); timerID = tic;
    for j=1:num_bones
        fprintf(', %d', j);
        bone_mask = ismember(mask_img,bone_label_indx{j});
        if(~any(bone_mask(:))), continue; end   % no target bone
        masked_volume = zeros(size(myu_img),'like',myu_img);
        masked_volume(bone_mask) = myu_img(bone_mask);

        if(mask_filter_width>0)
          masked_volume = ApplyEdgeBlurredMask(masked_volume, bone_mask, mask_filter_width);
%             masked_volume = ApplyErodedMask(masked_volume, bone_mask, mask_filter_width);
        end
        bb = GetBoundingBox(masked_volume,0);
        s = size(masked_volume);
        volume_offset_4x4xN(:,:,j) = [eye(3) -([s(1)/2-(bb(1)+bb(2)-1)/2, s(2)/2-(bb(3)+bb(4)-1)/2, s(3)/2-(bb(5)+bb(6)-1)/2] .* ElementSpacing)'; 0 0 0 1];

        masked_cropped_volume = masked_volume(bb(1):bb(2),bb(3):bb(4),bb(5):bb(6));
        volumePlan = regTools.CreateInterpolatorPlan( masked_cropped_volume, ElementSpacing );
        for i=1:size(DS,1)
            all_ElementSpacing_volume(i,:) = ElementSpacing.*DS(i,1); 
            regTools.SetVolumeInfo( struct('VolumeDim', ceil(size(masked_cropped_volume)./DS(i,1)), 'VoxelSize', all_ElementSpacing_volume(i,:)) );
            all_masked_volumes{i,j} = regTools.Interpolation( volumePlan, [0 0 0 0 0 0], RegTools.InterpolatorType_Bicubic, 0, -0.5 ); % bi-cubic interpolation
        end
        regTools.DeleteInterpolatorPlan( volumePlan );
    end
    fprintf(', done in %f sec\n', toc(timerID));
    