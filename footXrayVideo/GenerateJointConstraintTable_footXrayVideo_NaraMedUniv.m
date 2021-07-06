function joint_constraint_tbl = GenerateJointConstraintTable_footXrayVideo_NaraMedUniv( bone_mode, volume_offset_4x4xN, joint_center_tbl, bone_tbl )

switch bone_mode
    case {'tarsal_bones_proximal', 'all_tarsal_bones'}
        if(~all(ismember({'tibia','talus'},bone_tbl.Row)))
            joint_constraint_tbl = [];
            return
        end

        base_object = [bone_tbl{'talus',1}]';
        ref_object = [bone_tbl{'tibia',1}]';
        joint_center_label = {'tibia'};
        constraint_threshold = 5;
        constraint_type = repmat({'spherical'}, length(base_object), 1);
        
    case {'tibia', 'talus', 'calcaneus', 'navicular'}
        base_object = [];
        
    case {'tarsal_bones_distal'}
        base_object = [];
        
    case {'metatarsal_bones'}
        base_object = [];
        
    case 'all_metatarsal_bones'
        if(~all(ismember({'Medial_cuneiform','Intermediate_cuneiform','Lateral_cuneiform','cuboid','1st_metatarsal','2nd_metatarsal','3rd_metatarsal','4th_metatarsal','5th_metatarsal'},bone_tbl.Row)))
            joint_constraint_tbl = [];
            return
        end

        base_object = [bone_tbl{'talus',1}, bone_tbl{'Medial_cuneiform',1}, bone_tbl{'Intermediate_cuneiform',1}, bone_tbl{'Lateral_cuneiform',1}, bone_tbl{'cuboid',1}, bone_tbl{'cuboid',1}]';
        ref_object = [bone_tbl{'tibia',1}, bone_tbl{'1st_metatarsal',1}, bone_tbl{'2nd_metatarsal',1}, bone_tbl{'3rd_metatarsal',1}, bone_tbl{'4th_metatarsal',1}, bone_tbl{'5th_metatarsal',1}]';
        joint_center_label = {'tibia', '1st_metatarsal', '2nd_metatarsal', '3rd_metatarsal', '4th_metatarsal' '5th_metatarsal'};
        constraint_threshold = [5; zeros(length(base_object)-1,1)]; %ones(length(base_object),1) * 3;
        constraint_type = repmat({'spherical'}, length(base_object), 1);
end

if(~isempty(base_object))
    coordinate_4x4_wrt_base = cell(length(base_object),1);
    coordinate_4x4_wrt_ref = cell(length(ref_object),1);
    for i=1:length(base_object)
        coordinate_4x4_wrt_base{i} = volume_offset_4x4xN(:,:,base_object(i)) \ [eye(3) joint_center_tbl{joint_center_label{i},:}(:); 0 0 0 1];
        coordinate_4x4_wrt_ref{i} = volume_offset_4x4xN(:,:,ref_object(i)) \ [eye(3) joint_center_tbl{joint_center_label{i},:}(:); 0 0 0 1];
    end

    joint_constraint_tbl = table(constraint_type, base_object, ref_object, coordinate_4x4_wrt_base, coordinate_4x4_wrt_ref, constraint_threshold);
else
    joint_constraint_tbl = [];
end
