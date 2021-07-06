function img_2D_original = ApplyLCN_image_wise(img_2D_original, LCN_sigma)
if(LCN_sigma>0)
    for i=1:size(img_2D_original,4)
        for j=1:size(img_2D_original,3)
            img_2D_original(:,:,j,i) = localnormalize(img_2D_original(:,:,j,i), LCN_sigma, LCN_sigma);
        end
    end
    img_2D_original(isnan(img_2D_original)) = 0;
end
