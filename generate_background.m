function generate_background(f,p,ZB,Lx,Ly,dot_size_pixels,dot_spacing)

    dot_diameter = p*ZB*dot_size_pixels/f; 

    nx = Lx / ((dot_spacing+1)*dot_diameter);
    ny = Ly / ((dot_spacing+1)*dot_diameter);
    n = floor(nx*ny);
    dot_coordinates = [Lx*rand(1,n); Ly*rand(1,n)];
    
    figure('Units', 'pixels', 'Position', [100, 100, Lx*500/Ly, 500]); 
    axes('Position', [0, 0, 1, 1], 'Units', 'normalized');
    hold on; axis equal
    xlim([0, Lx]); ylim([0, Ly]);
    set(gca, 'XColor', 'none', 'YColor', 'none', 'Color', 'w')
    drawnow; 
    
    axis_pixels = gca().Position(3:4) .* gcf().Position(3:4);
    dot_diameter_pixels_x = dot_diameter * (axis_pixels(1) / Lx);
    dot_diameter_pixels_y = dot_diameter * (axis_pixels(2) / Ly);
    dot_diameter_pixels = mean([dot_diameter_pixels_x, dot_diameter_pixels_y]);
    
    screenDPI = get(0, 'ScreenPixelsPerInch');
    d_points = dot_diameter_pixels * (72 / screenDPI); 
    size =  d_points^2; 
    
    scatter(dot_coordinates(1,:),dot_coordinates(2,:),size,"k","filled");
    exportgraphics(gcf, 'background.pdf','ContentType','vector','Resolution',600)
    
    close

end

