video1_filename = '\\scallop\User\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data_model\image\Xray_interpolated\MODEL1\lateral\BP_LL0038_interpolated.mhd';
video2_filename = '\\scallop\User\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data_model\image\Xray_interpolated\MODEL1\oblique\BP_LL0018_interpolated.mhd';
output_filename = 'bone_model_Xray_interpolated';

% video1_filename = '\\scallop\User\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data_model\image\Xray_interpolated\MODEL1\lateral\BP_LL0038_inpainted_interpolated.mhd';
% video2_filename = '\\scallop\User\Nara Medic Univ Orthopaedic\Projects\Foot2D3D\data_model\image\Xray_interpolated\MODEL1\oblique\BP_LL0018_inpainted_interpolated.mhd';
% output_filename = 'bone_model_Xray_inpainted_interpolated';

output_dir = 'D:\Collaboration\Nara Medical University\20210131_Foot2D3D_summary';
if(~exist(output_dir, 'dir')), mkdir(output_dir); end

[video1, hdr1] = mhdread(video1_filename);
[video2, hdr2] = mhdread(video2_filename);

%%
figure('Position',[100 150 [1600 800]*3/10], 'PaperPositionMode', 'auto', 'Color', 'w'); colormap(gray(256));
writerObj = VideoWriter( fullfile(output_dir, output_filename) );
writerObj.FrameRate = 3;
open(writerObj);

for i=1:5:size(video2,3)
    subplot('Position',[0 0 0.5 1]);
    im(video1(:,:,i)); axis off; set(gca,'clim',[100 140]);
    subplot('Position',[0.5 0 0.5 1]);
    im(video2(:,:,i)); axis off; set(gca,'clim',[80 120]);
    writeVideo(writerObj, getframe(gcf));
end
writerObj = [];
