function floating = GetMergedFloatingImages( regTools, volumePlans, trans_param_4x4, num_views, target_object )

if(exist('target_object','var') && ~isempty(target_object))
%     trans_param_4x4 = trans_param_4x4(:,:,target_object,:,:);
    volumePlans = volumePlans(target_object);
end

for i=1:size(trans_param_4x4,3)
    if(i==1), MemoryStoreMode = regTools.MemoryStoreMode_Replace; else, MemoryStoreMode = regTools.MemoryStoreMode_Additive; end
    if(nargout==0)
        regTools.ForwardProject(volumePlans(i), reshape(trans_param_4x4(:,:,i,:,:),4,4,[]), [], num_views, [], [], MemoryStoreMode);
    else
        if(i<size(trans_param_4x4,3))
            regTools.ForwardProject(volumePlans(i), reshape(trans_param_4x4(:,:,i,:,:),4,4,[]), [], num_views, [], [], MemoryStoreMode);
        else
            floating = regTools.ForwardProject(volumePlans(i), reshape(trans_param_4x4(:,:,i,:,:),4,4,[]), [], num_views, [], [], MemoryStoreMode);
        end
    end
end
