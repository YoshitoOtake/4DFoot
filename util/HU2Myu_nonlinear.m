function myu_volume = HU2Myu_nonlinear(HU_volume, myu_water, HU_filter_threshold)

if(~exist('HU_filter_threshold','var')), HU_filter_threshold = 100; end
% myu_volume = max( ((1000+single(HU_volume))/2000).^3 * myu_water*2, 0 ); % crosses with linear model at +1000HU
myu_volume = max( ((max(single(HU_volume),0))/1000).^2 * myu_water*2, 0 ); % crosses with linear model at +1000HU
% myu_volume = max( ((max(single(HU_volume),0))/1000) * myu_water*2, 0 ); % crosses with linear model at +1000HU
% myu_volume = max( ((1000+single(HU_volume))/1000) * myu_water, 0 ); % linear model
myu_volume(HU_volume<HU_filter_threshold) = 0;

if(0)
    %%
    % sample test
    HU = -1000:1500;
    myu1 = max( ((1000+single(HU))/1000)    * myu_water  , 0 );
    myu2 = max( ((max(single(HU),0))/1000).^2 * myu_water*2, 0 );
    myu3 = max( ((1000+single(HU))/1000)    * myu_water  , 0 ); myu3(HU<0) = 0;
    figure; plot(HU, [myu1; myu2; myu3]);
end