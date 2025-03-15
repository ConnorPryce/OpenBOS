function u = shift_solution(u,u0,i,j)
    
    u = u - u(i,j,:) + u0;
 
end

