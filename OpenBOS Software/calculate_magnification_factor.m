function M = calculate_magnification_factor(Iref,Lref)
    
    figure('Name', 'OpenBOS - Length Calibration', ...
           'NumberTitle', 'off', ...
           'Position', [600, 200, 700, 500]);

    axes('Parent', gcf, ...
         'Position', [0.1, 0.3, 0.8, 0.65]);

    imshow(Iref)
    roi = drawline("Color","r");
    points = roi.Position;
    close
    Lpx = sqrt((points(2,1)-points(1,1))^2+(points(2,2)-points(1,2))^2);

    M = Lref/Lpx;

end

