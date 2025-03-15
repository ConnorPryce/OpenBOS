function mask = create_mask(frames)

    [ny,nx,n_frames] = size(frames);
    current_frame = 1;
    poly_handle = [];
    mask = ones(ny,nx,"logical");

    figure_handle = figure('Name', 'OpenBOS - Create Mask', ...
                           'NumberTitle', 'off', ...
                           'Position', [600, 200, 700, 500]);

    axis_handle = axes('Parent', figure_handle, ...
                       'Position', [0.1, 0.3, 0.8, 0.65]);
    
    image_handle = imshow(frames(:,:,current_frame),'Parent',axis_handle);
    
    slider_handle = uicontrol('Style', 'slider', ...
                              'Min', 1, ...
                              'Max', n_frames, ...
                              'Value', 1, ...
                              'SliderStep', [1/(n_frames-1), 1/(n_frames-1)], ...
                              'Units', 'normalized', ...
                              'Position', [0.15, 0.2, 0.7, 0.05]);

    count_handle = uicontrol('Style', 'text', ...
                             'String', 'Frame: 1', ...
                             'Units', 'normalized', ...
                             'Position', [0.15, 0.095, 0.15, 0.05], ...
                             'FontSize', 12,...
                             'HorizontalAlignment', 'left');

    mask_handle = uicontrol('Style', 'pushbutton', ...
                            'String', 'Mask', ...
                            'Units', 'normalized', ...
                            'Position', [0.425, 0.1, 0.15, 0.05], ...
                            'FontSize', 12, ...
                            'Callback',@mask_callback);

    save_handle = uicontrol('Style', 'pushbutton', ...
                            'String', 'Save', ...
                            'Units', 'normalized', ...
                            'Position', [0.7, 0.1, 0.15, 0.05], ...
                            'FontSize', 12, ...
                            'Callback',@save_callback);

    slider_handle.Callback =  @slider_callback;
    save_handle.Callback =  @save_callback;

    function update()
        count_handle.String = sprintf('Frame: %d', current_frame);
        image_handle.CData = frames(:,:,current_frame);
    end
    
    function slider_callback(~,~)
        current_frame = round(slider_handle.Value);
        update;
    end 
    
    function mask_callback(~,~)
        poly_handle = drawpolygon(axis_handle,'color','r');
    end 

    function save_callback(~,~)
        if isempty(poly_handle)
            mask = ones(nx,ny);
        else
            mask = ~createMask(poly_handle);
            mask = image_to_ndgrid(mask);
        end 
           
        close
    end 

    waitfor(figure_handle);

end

