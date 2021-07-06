function masked_volume = ApplyEdgeBlurredMask(volume, mask, blur_width)

is2D = 0;

if(blur_width == 0)
    masked_volume = volume .* single(mask);
    return;
end

% tic;
if(is2D)
    normalized_dist2D = zeros(size(mask), 'single');
    non_zero_slice = find(squeeze(sum(sum(mask,1),2))>0);
    for i=1:length(non_zero_slice)
        normalized_dist2D(:,:,non_zero_slice(i)) = 1-mat2gray(bwdist(mask(:,:,non_zero_slice(i))),[0 blur_width]);
    end
    masked_volume = volume.*normalized_dist2D;
else
    bb = GetBoundingBox(mask, 0) + [-1 1 -1 1 -1 1]*blur_width;
    bb([1 3 5]) = max(bb([1 3 5]),1);
    bb([2 4 6]) = min(bb([2 4 6]),size(mask));
    mask_erode = imerode(mask(bb(1):bb(2),bb(3):bb(4),bb(5):bb(6)),strel('cube',blur_width));
    normalized_dist3D = 1-mat2gray(bwdist(mask_erode), [0 blur_width]);
    masked_volume = zeros(size(mask),'single');
    masked_volume(bb(1):bb(2),bb(3):bb(4),bb(5):bb(6)) = volume(bb(1):bb(2),bb(3):bb(4),bb(5):bb(6)).*normalized_dist3D;
end
% fprintf('done in %f sec on CPU\n', toc);
