function pos2D_Mx2xN = cleanup_2D_landmark(pos2D_Mx2xN)

[M, ~, N] = size(pos2D_Mx2xN);  % M: number of landmarks, N: number of frames
x = 1:N;
for i=1:M
    for j=1:2
        % Obtain NaN index
        nanidx = ~isnan(squeeze(pos2D_Mx2xN(i,j,:)));
        if(~isempty(nanidx))
            % inter(extra-)polate NaN data
            pos2D_Mx2xN(i,j,~nanidx) = interp1(x(nanidx), squeeze(pos2D_Mx2xN(i,j,nanidx)), x(~nanidx), 'pchip');
        end
    end
end