function joint_constraint_tbl = GenerateJointConstraintTable_subtalar_joint( volume_offset_4x4xN, joint_center_tbl, bone_tbl )

if(~all(ismember({'tibia','talus'},bone_tbl.Row)))
    joint_constraint_tbl = [];
    return
end

base_object = [bone_tbl{'talus',1}]';
ref_object = [bone_tbl{'tibia',1}]';
joint_center_label = {'tibia'};
constraint_threshold = [5]';
constraint_type = repmat({'spherical'}, length(base_object), 1);

coordinate_4x4_wrt_base = cell(length(base_object),1);
coordinate_4x4_wrt_ref = cell(length(ref_object),1);
for i=1:length(base_object)
    coordinate_4x4_wrt_base{i} = volume_offset_4x4xN(:,:,base_object(i)) \ [eye(3) joint_center_tbl{joint_center_label{i},:}(:); 0 0 0 1];
    coordinate_4x4_wrt_ref{i} = volume_offset_4x4xN(:,:,ref_object(i)) \ [eye(3) joint_center_tbl{joint_center_label{i},:}(:); 0 0 0 1];
end
joint_constraint_tbl = table(constraint_type, base_object, ref_object, coordinate_4x4_wrt_base, coordinate_4x4_wrt_ref, constraint_threshold);
