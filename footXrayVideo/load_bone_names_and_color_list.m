% object_names = {'tibia', 'fibula', 'talus', 'calc', 'navi'};
% object_colors = repmat([241 214 145]/255, [length(object_names),1]);

object_names = {
    'tibia','fibula','talus','calc','navi','Medial_cuneiform','Intermediate_cuneiform','Lateral_cuneiform','cuboid', ...
    '1st_metatarsal','2rd_metatarsal','3rd_metatarsal','4th_metatarsal','5th_metatarsal', ...
    '1st_phalanx_proximalis','2rd_phalanx_proximalis','3rd_phalanx_proximalis','4th_phalanx_proximalis','5th_phalanx_proximalis'
};

object_colors = [
    [241 214 145]/255-[0 0.3 0.3]; [241 214 145]/255; 0 1 1; 0.75 1 0.25;
    1 1 0; 0 1 0; 1 0.5 0.5; 1 0.5 0.5; 0.5 0 0.5;
    0 0 1; 1 0 0; 1 0 1; 1 0.5 0; 0 1 1;
    1 0 0; 1 1 0; 1 0.5 0; 1 0 1; 0 0 1; ];
