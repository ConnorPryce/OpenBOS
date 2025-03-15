function [dndx,dndy] = refraction_gradients(dx,dy,n0,f,ZA,ZD,D,p)
    
    dndx = n0*(ZA+ZD)*dx*p/(ZD*f*D);
    dndy = n0*(ZA+ZD)*dy*p/(ZD*f*D);

end

