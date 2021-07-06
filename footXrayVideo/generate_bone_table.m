function bone_tbl = generate_bone_table( bone_mode )

bone_list = jsondecode(fileread('ankle_new.json'));
bone_names = cellfun(@(x) x{2},bone_list.color,'UniformOutput',false);

switch bone_mode
    case 'tarsal_bones_proximal'
        % do not include metatarsal bones
        bone_label = bone_names([1 3:5]); %{'tibia', 'talus', 'calcaneal', 'navicular''}';
        bone_label_indx = {[1,2], 3, 4, 5}';
    case 'tarsal_bones_distal'
        % do not include metatarsal bones
        bone_label = bone_names([6:9]); %{'Medial_cuneiform', 'Intermediate_cuneiform', 'Lateral_cuneiform', 'cuboid' '}';
        bone_label_indx = {6, 7, 8, 9}';
    case 'metatarsal_bones'
        % include metatarsal bones
        bone_label = bone_names([10:14]); %{'tibia', 'talus', 'calcaneal', 'navicular', 'Medial_cuneiform', 'Intermediate_cuneiform', 'Lateral_cuneiform', 'cuboid', '1st_metatarsal', '}';
        bone_label{contains(bone_label,'2rd_metatarsal')} = '2nd_metatarsal';
        bone_label_indx = {10, 11, 12, 13, 14}';
    case 'tibia'
        bone_label = bone_names(1);        bone_label_indx = {[1,2]}';
    case 'talus'
        bone_label = bone_names(3);        bone_label_indx = {3}';
    case 'calcaneus'
        bone_label = bone_names(4);        bone_label_indx = {4}';
    case 'navicular'
        bone_label = bone_names(5);        bone_label_indx = {5}';
    case 'all_tarsal_bones'
        % include metatarsal bones
%         bone_label = bone_names([1 3:9]); %{'tibia', 'talus', 'calcaneal', 'navicular', 'Medial_cuneiform', 'Intermediate_cuneiform', 'Lateral_cuneiform', 'cuboid' '}';
%         bone_label_indx = {[1,2], 3, 4, 5, 6, 7, 8, 9}';
        bone_label = bone_names([1 3:5 9]); %{'tibia', 'talus', 'calcaneal', 'navicular', 'Medial_cuneiform', 'Intermediate_cuneiform', 'Lateral_cuneiform', 'cuboid' '}';
        bone_label_indx = {[1,2], 3, 4, 5, 9}';
    case 'all_metatarsal_bones'
        % include metatarsal bones
        bone_label = bone_names([1 3:14]); %{'tibia', 'talus', 'calcaneal', 'navicular', 'Medial_cuneiform', 'Intermediate_cuneiform', 'Lateral_cuneiform', 'cuboid', '1st_metatarsal', '}';
        bone_label{contains(bone_label,'2rd_metatarsal')} = '2nd_metatarsal';
        bone_label_indx = {[1,2], 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14}';
end

bone_tbl = table((1:length(bone_label))', bone_label_indx, 'VariableNames', {'indx', 'label_indx'}, 'RowNames', bone_label);    % include metatarsal bones
