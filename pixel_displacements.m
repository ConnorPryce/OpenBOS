function [X,Y,ds,gx,gy,domain_mask] = pixel_displacements(Ir,I,window_size,overlap,image_mask)
    
    Ir = single(image_to_ndgrid(Ir));
    I = single(image_to_ndgrid(I));

    [npx,npy,nframes] = size(I);
    ds = window_size-overlap;
    nwx = floor((npx-overlap)/ds);
    nwy = floor((npy-overlap)/ds);
    half_window_size = floor(window_size/2) + 1;
    x = single([half_window_size, half_window_size + (1:nwx-1)*ds]);
    y = single([half_window_size, half_window_size + (1:nwy-1)*ds]);
    [X,Y] = ndgrid(x,y);
   
    valid_region = conv2(image_mask,ones(max([window_size,5])),'same') == window_size^2;
    domain_mask = valid_region(x,y);

    sigma = window_size/6;
    weight = fspecial('gaussian', window_size, sigma); 

    kernal = single([1/12,-2/3,0,2/3,-1/12])';
    
    Irx = conv2(Ir, kernal, 'same');
    Iry = conv2(Ir, kernal', 'same');
    Irt = conv2(Ir, -ones(3,3)/9, 'same');

    gx = nan(nwx,nwy,nframes,"single");
    gy = nan(nwx,nwy,nframes,"single");

    parfor frame = 1:nframes
        
        Ix = 0.5 * (Irx + conv2(I(:,:,frame), kernal, 'same'));
        Iy = 0.5 * (Iry + conv2(I(:,:,frame), kernal', 'same'));
        It = (Irt + conv2(I(:,:,frame), ones(3,3)/9, 'same'));

        IxIxSum = conv2(Ix.^2, weight, 'same');
        IyIySum = conv2(Iy.^2, weight, 'same');
        IxIySum = conv2(Ix.*Iy, weight, 'same');
        IxItSum = conv2(Ix.*It, weight, 'same');
        IyItSum = conv2(Iy.*It, weight, 'same');

        temp_dx = nan(nwx,nwy,"single");
        temp_dy = nan(nwx,nwy,"single");
    
        for i = 1:nwx
            for j = 1:nwy
                ix = x(i);
                iy = y(j);
                if domain_mask(i,j)
                    D = IxIxSum(ix, iy) * IyIySum(ix, iy) - IxIySum(ix, iy)^2;                 
                    if D > 1e-3
                        temp_dx(i, j) = (-IxItSum(ix, iy) * IyIySum(ix, iy) + IyItSum(ix, iy) * IxIySum(ix, iy)) / D;
                        temp_dy(i, j) = (-IxIxSum(ix, iy) * IyItSum(ix, iy) + IxIySum(ix, iy) * IxItSum(ix, iy)) / D;
                    else
                        temp_dx(i, j) = 0;
                        temp_dy(i, j) = 0;
                    end
                end 
            end
        end
        
        gx(:, :, frame) = temp_dx;
        gy(:, :, frame) = temp_dy;
        
    end
end
