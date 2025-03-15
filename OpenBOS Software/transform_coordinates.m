function [Xt,Yt,dst] = transform_coordinates(M,X,Y,ds)
    
    Xt = M*X;
    Yt = M*Y;
    dst = M*ds;
    
end

