function out_fv = merge_fv(in1_fv, in2_fv)

offset = size(in1_fv.vertices,1);
out_fv = struct('faces', [in1_fv.faces; in2_fv.faces + offset], 'vertices', [in1_fv.vertices; in2_fv.vertices]);