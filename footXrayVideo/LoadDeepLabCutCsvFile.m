function [positions, indx, likelihood] = LoadDeepLabCutCsvFile( csv_filename )

fid = fopen(csv_filename, 'r');
fgetl(fid);
fgetl(fid);
text = textscan(fgetl(fid),'%s','delimiter',',');
x_indx = strcmp(text{1},'x');
y_indx = strcmp(text{1},'y');
likelihood_indx = strcmp(text{1},'likelihood');
positions = [];
likelihood = [];
indx = {};
while(~feof(fid))
    text = textscan(fgetl(fid),'%s','delimiter',',');
    if(length(text{1})<1), continue; end    % skip empty line
    if(length(text{1})<length(x_indx)), text{1}(length(x_indx)) = {[]}; end
    positions_x = str2double(text{1}(x_indx));
    positions_y = str2double(text{1}(y_indx));
    positions = cat(3, positions, [positions_x positions_y]);
    likelihood_one_line = str2double(text{1}(likelihood_indx));
    likelihood = cat(2, likelihood, likelihood_one_line);
    indx = cat(1, indx, text{1}(1));
end
fclose(fid);

% tbl = readtable(csv_filename,'ReadVariableNames',false,'VariableNamingRule','preserve');
% positions_x = str2double(tbl{4:end,strcmp(tbl{3,:},'x')})';
% positions_y = str2double(tbl{4:end,strcmp(tbl{3,:},'y')})';
% positions = cat(2, reshape(positions_x,[],1,height(tbl)-3), reshape(positions_y,[],1,height(tbl)-3));
% likelihood = str2double(tbl{4:end,strcmp(tbl{3,:},'likelihood')});
% 
% indx = tbl{4:end,1};
