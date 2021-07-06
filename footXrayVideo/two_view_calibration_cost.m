function cost = two_view_calibration_cost(param, opt_data, visualization_on, optimization_mode, force_visualize)
if(~exist('force_visualize','var')), force_visualize = false; end
    
    M = size(param,2);
    cost = zeros(1, M);
    
%     if(optimization_mode==4),    disp(param(1,:)); end
    scaled_param = opt_data.initial_param(:)*ones(1,M) + diag(opt_data.param_scaling)*param;
    switch optimization_mode
        case 1
            % for 1-view, 1-image -> 6-DOF optimization
            ProjectionMatrix1 = [[-opt_data.focal_length_pix(1), 0, opt_data.origin_2D_pix(1); 0, -opt_data.focal_length_pix(1), opt_data.origin_2D_pix(2); 0, 0, 1] [0 0 0]'];
            num_box_positions = 1;
            all_dist_errors = cell(num_box_positions, 1);
        case 2
            % for 1-view, 6-DOF box position + focal length -> 7-DOF
            % optimization (seems like unstable for interaction between
            % focal length and z-position)
            num_box_positions = 1;
            all_dist_errors = cell(num_box_positions, 1);
        case 3
            % for 2-view, optimize only relative camera position -> 6-DOF optimization
            ProjectionMatrix1 = [[-opt_data.focal_length_pix(1), 0, opt_data.origin_2D_pix(1); 0, -opt_data.focal_length_pix(1), opt_data.origin_2D_pix(2); 0, 0, 1] [0 0 0]'];
            ProjectionMatrix2 = [[-opt_data.focal_length_pix(2), 0, opt_data.origin_2D_pix(1); 0, -opt_data.focal_length_pix(2), opt_data.origin_2D_pix(2); 0, 0, 1] [0 0 0]'];
            PM_offset = 0;
            num_box_positions = size(opt_data.calib_box_positions,2); %(size(param,1)-6)/6;
            all_dist_errors = cell(num_box_positions, 2);
        case 4
            % for 2-view, optimize relative camera position + focal length + origin of each view -> 6-DOF optimization
            PM_offset = 6;
            num_box_positions = size(opt_data.calib_box_positions,2); %(size(param,1)-6)/6;
            all_dist_errors = cell(num_box_positions, 2);
    end
    for i=1:M
        if(optimization_mode==4)
            ProjectionMatrix1 = [[-scaled_param(1,i), 0, scaled_param(3,i); 0, -scaled_param(1,i), scaled_param(4,i); 0, 0, 1] [0 0 0]'];
            ProjectionMatrix2 = [[-scaled_param(2,i), 0, scaled_param(5,i); 0, -scaled_param(2,i), scaled_param(6,i); 0, 0, 1] [0 0 0]'];
        end
        for j=1:num_box_positions
            switch optimization_mode
                case 1
                    T = RegTools.convertTransRotTo4x4( scaled_param(1:6,i) );
                case 2
                    T = RegTools.convertTransRotTo4x4( scaled_param(1+(1:6),i) );
                    ProjectionMatrix1 = [[-scaled_param(1,i), 0, opt_data.origin_2D_pix(1); 0, -scaled_param(1,i), opt_data.origin_2D_pix(2); 0, 0, 1] [0 0 0]'];
                case 3
                    T = RegTools.convertTransRotTo4x4( opt_data.calib_box_positions(:,j) ); %RegTools.convertTransRotTo4x4( scaled_param((PM_offset+6+(j-1)*6)+(1:6),i) );
                    Camera2pos = opt_data.CameraPosOffset(:,:,1) * RegTools.convertTransRotTo4x4( scaled_param(PM_offset + (1:6),i) ) * opt_data.CameraPosOffset(:,:,2);
                    p2 = ProjectionMatrix2 * inv(Camera2pos) * T * [opt_data.pos_3D'; ones(1,size(opt_data.pos_3D,1))];
                    projected2 = (p2(1:2,:)./[p2(3,:); p2(3,:)])';
                    all_dist_errors{j,2} = sqrt(sum((projected2-opt_data.pos_2D{j,2}).^2,2));
                case 4
%                     T = RegTools.convertTransRotTo4x4( opt_data.calib_box_positions(:,j) ); %RegTools.convertTransRotTo4x4( scaled_param((PM_offset+6+(j-1)*6)+(1:6),i) );
                    T = RegTools.convertTransRotTo4x4( scaled_param((PM_offset+6+(j-1)*6)+(1:6),i) );
                    Camera2pos = opt_data.CameraPosOffset(:,:,1) * RegTools.convertTransRotTo4x4( scaled_param(PM_offset + (1:6),i) ) * opt_data.CameraPosOffset(:,:,2);
                    p2 = ProjectionMatrix2 * inv(Camera2pos) * T * [opt_data.pos_3D'; ones(1,size(opt_data.pos_3D,1))];
                    projected2 = (p2(1:2,:)./[p2(3,:); p2(3,:)])';
                    all_dist_errors{j,2} = sqrt(sum((projected2-opt_data.pos_2D{j,2}).^2,2));
            end
            p1 = ProjectionMatrix1 * T * [opt_data.pos_3D'; ones(1,size(opt_data.pos_3D,1))];
            projected1 = (p1(1:2,:)./[p1(3,:); p1(3,:)])';
            all_dist_errors{j,1} = sqrt(sum((projected1-opt_data.pos_2D{j,1}).^2,2));
        end
        cost(i) = mean(cell2mat(reshape(all_dist_errors,[],1)));    % mean reprojection error
    end
    opt_data.cost_log = [opt_data.cost_log; cost'];
    opt_data.RegistrationLogs = cat(1, opt_data.RegistrationLogs, scaled_param');


    if(visualization_on)
        if(opt_data.iteration_count == 0 || opt_data.iteration_count>(opt_data.previous_show_iteration + opt_data.progress_show_step) || force_visualize)
            [min_cost, min_cost_indx] = min(cost, [], 2);
            scaled_param = opt_data.initial_param(:) + param(:,min_cost_indx) .* opt_data.param_scaling(:);
            switch optimization_mode
                case 1
                    focal_length_mm = opt_data.focal_length_pix * mean(opt_data.img2D_ElementSpacing);
                    image_plane_centers = opt_data.origin_2D_pix;
                    all_Ts = zeros(4, 4, 1);
                case 2
                    focal_length_mm = scaled_param(1) * mean(opt_data.img2D_ElementSpacing);
                    image_plane_centers = opt_data.origin_2D_pix;
                    all_Ts = zeros(4, 4, 1);
                case 3
                    focal_length_mm = [opt_data.focal_length_pix(1)*mean(opt_data.img2D_ElementSpacing) opt_data.focal_length_pix(2)*mean(opt_data.img2D_ElementSpacing)];
                    image_plane_centers = [opt_data.origin_2D_pix; opt_data.origin_2D_pix];
                    all_Ts = zeros(4, 4, num_box_positions);
                case 4
                    focal_length_mm = [scaled_param(1)*mean(opt_data.img2D_ElementSpacing) scaled_param(2)*mean(opt_data.img2D_ElementSpacing)];
                    image_plane_centers = [scaled_param(3:4); scaled_param(5:6)];
                    all_Ts = zeros(4, 4, num_box_positions);
            end
            image_plane_corners = opt_data.origin_2D_pix([1 1 2 2])*2 + [0 1 0 1];          
            
            moved_3D = zeros(8, 3, num_box_positions);
            projected = cell(num_box_positions, 2);
            opt_data.landmark_errors = cell(num_box_positions, 2);
            for j=1:num_box_positions
                switch optimization_mode
                    case 1
                        all_Ts(:,:,j) = RegTools.convertTransRotTo4x4( scaled_param(1:6) );
                    case 2
                        all_Ts(:,:,j) = RegTools.convertTransRotTo4x4( scaled_param(1+(1:6)) );
                        ProjectionMatrix1 = [[-scaled_param(1), 0, opt_data.origin_2D_pix(1); 0, -scaled_param(1), opt_data.origin_2D_pix(2); 0, 0, 1] [0 0 0]'];
                    case 3
                        all_Ts(:,:,j) = RegTools.convertTransRotTo4x4( opt_data.calib_box_positions(:,j) ); %RegTools.convertTransRotTo4x4( scaled_param((PM_offset+6+(j-1)*6)+(1:6)) );
                        Camera2pos = opt_data.CameraPosOffset(:,:,1) * RegTools.convertTransRotTo4x4( scaled_param(PM_offset + (1:6)) ) * opt_data.CameraPosOffset(:,:,2);
                        p2 = ProjectionMatrix2 * inv(Camera2pos) * all_Ts(:,:,j) * [opt_data.pos_3D'; ones(1,size(opt_data.pos_3D,1))];
                        projected{j,2} = (p2(1:2,:)./[p2(3,:); p2(3,:)])';
                        opt_data.landmark_errors{j,2} = sqrt(sum((projected{j,2}-opt_data.pos_2D{j,2}).^2,2));
                    case 4
                        ProjectionMatrix1 = [[-scaled_param(1), 0, scaled_param(3); 0, -scaled_param(1), scaled_param(4); 0, 0, 1] [0 0 0]'];
                        ProjectionMatrix2 = [[-scaled_param(2), 0, scaled_param(5); 0, -scaled_param(2), scaled_param(6); 0, 0, 1] [0 0 0]'];
                        
%                         all_Ts(:,:,j) = RegTools.convertTransRotTo4x4( opt_data.calib_box_positions(:,j) ); 
                        all_Ts(:,:,j) = RegTools.convertTransRotTo4x4( scaled_param((PM_offset+6+(j-1)*6)+(1:6)) );
                        Camera2pos = opt_data.CameraPosOffset(:,:,1) * RegTools.convertTransRotTo4x4( scaled_param(PM_offset + (1:6)) ) * opt_data.CameraPosOffset(:,:,2);
                        p2 = ProjectionMatrix2 * inv(Camera2pos) * all_Ts(:,:,j) * [opt_data.pos_3D'; ones(1,size(opt_data.pos_3D,1))];
                        projected{j,2} = (p2(1:2,:)./[p2(3,:); p2(3,:)])';
                        opt_data.landmark_errors{j,2} = sqrt(sum((projected{j,2}-opt_data.pos_2D{j,2}).^2,2));
               end                
                moved_3D(:,:,j) = ApplyTransform4x4(opt_data.pos_3D, all_Ts(:,:,j));
                p1 = ProjectionMatrix1 * all_Ts(:,:,j) * [opt_data.pos_3D'; ones(1,size(opt_data.pos_3D,1))];
                projected{j,1} = (p1(1:2,:)./[p1(3,:); p1(3,:)])';
                opt_data.landmark_errors{j,1} = sqrt(sum((projected{j,1}-opt_data.pos_2D{j,1}).^2,2));
            end
            
            % visualize current estimate
            [hs,vs,tb,bb,lb,rb] = deal(0.07, 0.07, 0.1, 0.1, 0.01, 0.01);
            clf; colormap(jet(256));
            
            msubplot(2,2,1,hs,vs,tb,bb,lb,rb);
            hold on;
            for j=1:num_box_positions
                scatter3(moved_3D(:,1,j),moved_3D(:,2,j),moved_3D(:,3,j),100,'r.');
                line(moved_3D(:,1,j), moved_3D(:,2,j), moved_3D(:,3,j), 'Color', 'r');
                draw_coordinate_system(all_Ts(:,:,j), 100, 3);
                if(j==1)
                    for i=1:size(moved_3D,1), text(moved_3D(i,1,j),moved_3D(i,2,j),moved_3D(i,3,j),sprintf('%d',i),'Color','r', 'FontSize', 14); end
                end
            end
            fv = cube_fv(opt_data.bounding_box_3D, all_Ts(:,:,1));
            patch(fv, 'FaceColor', 'none', 'EdgeColor', 'k');
            
            draw_coordinate_system(eye(4), 300, 3);
            [X,Y] = meshgrid(linspace(image_plane_corners(1),image_plane_corners(2),10),linspace(image_plane_corners(3),image_plane_corners(4),10));
            Z = ones(size(X))*(-focal_length_mm(1));
            mesh(X,Y,Z, 'FaceColor', 'none', 'EdgeColor', 'k');
            if(optimization_mode == 3)
                draw_coordinate_system(Camera2pos, 300, 3);
                Z = ones(size(X))*(-focal_length_mm(2));
                XYZ_moved = ApplyTransform4x4([X(:) Y(:) Z(:)], inv(Camera2pos));
                mesh(reshape(XYZ_moved(:,1),size(X)),reshape(XYZ_moved(:,2),size(X)),reshape(XYZ_moved(:,3),size(X)), 'FaceColor', 'none', 'EdgeColor', 'k');
            end
            axis equal;
            set(gca,'xlim',[-1000 1000],'ylim',[-500 500],'zlim',[-1500 100],'ydir','reverse','box','off');
            camtarget([0 0 -500]); campos([-700 200 -700]); camup([0 1 0]);
            set(gca,'xdir','reverse');
            xlabel('X axis'); ylabel('Y axis'); zlabel('Z axis');

            if(length(opt_data.cost_log)>1)
                msubplot(2,2,2,hs,vs,tb,bb,lb,rb);
                plot(opt_data.cost_log)
%                 ylim = [0 10];
%                 if(optimization_mode), ylim = [0 1]; end
                set(gca, 'FontSize', 16, 'xlim', [1 length(opt_data.cost_log)], 'YScale', 'log'); %, 'ylim', ylim);
                ylabel('mean reprojection error [mm]'); xlabel('# of iterations');
            end

%             colors = {'g', 'b', 'm'};
            colors = mat2cell(hsv(num_box_positions), ones(num_box_positions,1), 3);
            msubplot(2,2,3,hs,vs,tb,bb,lb,rb);
            hold on;
            for k=1:num_box_positions
                scatter(projected{k,1}(:,1),projected{k,1}(:,2),100,colors{k},'Marker','.');
                line(projected{k,1}(:,1), projected{k,1}(:,2), 'Color', colors{k});
                scatter(opt_data.pos_2D{k,1}(:,1),opt_data.pos_2D{k,1}(:,2),20,colors{k},'Marker','o');
                line(opt_data.pos_2D{k,1}(:,1), opt_data.pos_2D{k,1}(:,2), 'Color', colors{k}, 'LineStyle', '--', 'LineWidth', 2);
                for i=1:size(opt_data.pos_2D{k,1},1), text(opt_data.pos_2D{k,1}(i,1)+5,opt_data.pos_2D{k,1}(i,2),sprintf('%d',i),'Color',colors{k}, 'FontSize', 14); end
            end            
            set(gca,'xlim',[0 512],'ylim',[0 512],'ydir','reverse','box','on','DataAspectRatio',[1 1 1]);

            if(optimization_mode == 3 || optimization_mode == 4)
                msubplot(2,2,4,hs,vs,tb,bb,lb,rb);
                hold on;
                for k=1:num_box_positions
                    scatter(projected{k,2}(:,1),projected{k,2}(:,2),100,colors{k},'Marker','.');
                    line(projected{k,2}(:,1), projected{k,2}(:,2), 'Color', colors{k});
                    scatter(opt_data.pos_2D{k,2}(:,1),opt_data.pos_2D{k,2}(:,2),20,colors{k},'Marker','o');
                    line(opt_data.pos_2D{k,2}(:,1), opt_data.pos_2D{k,2}(:,2), 'Color', colors{k}, 'LineStyle', '--', 'LineWidth', 2);
                    for i=1:size(opt_data.pos_2D{k,2},1), text(opt_data.pos_2D{k,2}(i,1)+5,opt_data.pos_2D{k,2}(i,2),sprintf('%d',i),'Color',colors{k}, 'FontSize', 14); end
                end            
                set(gca,'xlim',[0 512],'ylim',[0 512],'ydir','reverse','box','on','DataAspectRatio',[1 1 1]);
            end
            
            btitle(sprintf('Iteration: %d, (%.3f, %.3f, %.3f,  %.3f, %.3f, %.3f), %d DOF, focal length: %f, cost: %f, %.3f sec', ...
                        opt_data.iteration_count, scaled_param(1:6), size(param,1), opt_data.focal_length_pix(1), min_cost, ...
                        toc(opt_data.TimerID)), 16, [0 0 0], 'none', 0.97);
            fprintf(sprintf('Iteration: %d, (%.3f, %.3f, %.3f,  %.3f, %.3f, %.3f), %d DOF, focal length: %f, cost: %f, %.3f sec\n', ...
                        opt_data.iteration_count, scaled_param(1:6), size(param,1), opt_data.focal_length_pix(1), min_cost, ...
                        toc(opt_data.TimerID)));
            drawnow;
            if(~isempty(opt_data.writerObj))
                writeVideo(opt_data.writerObj, getframe(gcf));
            end
            opt_data.previous_show_iteration = opt_data.iteration_count;
        end
        opt_data.iteration_count = opt_data.iteration_count + M;
    end
end

function r = rotx(t)
	ct = cos(t);
	st = sin(t);
	r =    [    1, 0, 0;   0, ct,-st;   0, st, ct];
end

function r = roty(t)
	ct = cos(t);
	st = sin(t);
	r =    [    ct	0	st;   0	1	0;   -st	0	ct];
end

function r = rotz(t)
	ct = cos(t);
	st = sin(t);
	r =    [ct	-st	0;  st	ct	0;  0	0	1];
end

function transformed = ApplyTransform4x4(points_Nx3, transform_4x4)
    temp = transform_4x4 * [points_Nx3'; ones(1,size(points_Nx3,1))];
    transformed = temp(1:3,:)';
end

function draw_coordinate_system( trans_4x4, axis_length, axis_width )
    vecs = trans_4x4(1:3,1:3) * axis_length;
    axis_color = {'r', 'g', 'b'};
    for i=1:3
        line(trans_4x4(1,4) + [0 vecs(1,i)], trans_4x4(2,4) + [0 vecs(2,i)], trans_4x4(3,4) + [0 vecs(3,i)], 'Color', axis_color{i}, 'LineWidth', axis_width);
    end
end

function fv = cube_fv(edges, transform_4x4)
    % cube_fv - return a patch structure (faces, vertices) of a 3D-cube
    %
    %   cube_fv(EDGES,ORIGIN,R) return a patch structure (faces, vertices) of a 3D-cube
    %   with the following properties:
    %   * EDGES : 3-element vector that defines the length of cube edges
    %   * transform_4x4: 4x4-element matrix that defines the transformation

    vertices = [0 0 0;1 0 0;1 1 0;0 1 0;0 0 1;1 0 1;1 1 1;0 1 1];
    faces = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
    fv = struct('faces', faces, 'vertices', ApplyTransform4x4(vertices*diag(edges), transform_4x4));
end
