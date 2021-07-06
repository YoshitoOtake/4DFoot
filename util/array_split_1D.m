function splitted = array_split_1D(array_1D, num_split)

num_output_elements = ceil(length(array_1D)/num_split);
residual = num_output_elements*num_split - length(array_1D);
num_output_elements_array = ones(num_split,1)*num_output_elements;
num_output_elements_array((end-residual+1):end) = num_output_elements_array((end-residual+1):end)-1;

splitted = mat2cell(array_1D(:), num_output_elements_array);
