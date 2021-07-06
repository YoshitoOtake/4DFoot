function draw_coordinate_system( trans_4x4, axis_length, axis_width )
    vecs = trans_4x4(1:3,1:3) * axis_length;
    axis_color = {'r', 'g', 'b'};
    for i=1:3
        line(trans_4x4(1,4) + [0 vecs(1,i)], trans_4x4(2,4) + [0 vecs(2,i)], trans_4x4(3,4) + [0 vecs(3,i)], 'Color', axis_color{i}, 'LineWidth', axis_width);
    end
end