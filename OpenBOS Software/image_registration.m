function [cx,cy,coeffs] = image_registration(X,Y,dx,dy,mask)
   
    nframes = size(dx,3);

    X = X - 0.5*(min(X(:)) + max(X(:)));
    Y = Y - 0.5*(min(Y(:)) + max(Y(:)));

    flat_mask = mask(:);
    X_masked = X(flat_mask);
    Y_masked = Y(flat_mask);
    
    n = sum(mask(:));
    Xsum = sum(X_masked);
    Ysum = sum(Y_masked);
    XYsum = sum(X_masked.*Y_masked);
    XXsum = sum(X_masked.*X_masked);
    YYsum = sum(Y_masked.*Y_masked);
    
    A = [   n,  Xsum,  Ysum;
         Xsum, XXsum, XYsum;
         Ysum, XYsum, YYsum];
    
    cx = zeros(size(dx));
    cy = zeros(size(dx));
    coeffs = zeros(6,nframes);

   parfor frame = 1:nframes

        Dx_masked = dx(:,:,frame);
        Dx_masked = Dx_masked(flat_mask);
        Dy_masked = dy(:,:,frame);
        Dy_masked = Dy_masked(flat_mask);
        
        Dxsum = sum(Dx_masked);
        Dysum = sum(Dy_masked);
        DxXsum = sum(Dx_masked.*X_masked);
        DxYsum = sum(Dx_masked.*Y_masked);
        DyXsum = sum(Dy_masked.*X_masked);
        DyYsum = sum(Dy_masked.*Y_masked);

        bx = [-Dxsum;-DxXsum;-DxYsum];
        by = [-Dysum;-DyXsum;-DyYsum];

        solx = A\bx;
        soly = A\by;
        
        cx(:,:,frame) = mask.*(solx(1) + solx(2)*X + solx(3)*Y);
        cy(:,:,frame) = mask.*(soly(1) + soly(2)*X + soly(3)*Y);
        coeffs(:,frame) = [solx;soly];
     
   end

end
