function masked_volume = ApplyErodedMask(volume, mask, erosion_width)
    bb = GetBoundingBox(mask, 0);
    bb([1 3 5]) = max(bb([1 3 5]),1);
    bb([2 4 6]) = min(bb([2 4 6]),size(mask));
    mask_erode = imerode(mask(bb(1):bb(2),bb(3):bb(4),bb(5):bb(6)),strel('cube',erosion_width));
    masked_volume = zeros(size(mask),'like',volume);
    masked_volume(bb(1):bb(2),bb(3):bb(4),bb(5):bb(6)) = volume(bb(1):bb(2),bb(3):bb(4),bb(5):bb(6)).*mask_erode;
