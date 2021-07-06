function file = getLatestTimestampedFile( query_wildcard )

local_reg_result_file = dir( query_wildcard );
if(isempty(local_reg_result_file))
    fprintf('no file (%s)\n', query_wildcard);
    file = [];
    return;
end
% no_ext_splitted = cellfun(@(x,y) regexp(x(1:y-1),'_','split'), {local_reg_result_file.name}',strfind({local_reg_result_file.name}', '.'), 'UniformOutput', false);
no_ext_splitted = cellfun(@(x,y) regexp(x,'_','split'), {local_reg_result_file.name}', 'UniformOutput', false);
timestamps = datetime(cellfun(@(x) [x{end-1} x{end}],no_ext_splitted,'UniformOutput',false),'InputFormat','yyyyMMddHHmmss');
[~,indx] = sort(timestamps);
file = local_reg_result_file(indx(end));

