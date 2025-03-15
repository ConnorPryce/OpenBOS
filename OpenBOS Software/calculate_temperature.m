function T = calculate_temperature(rho,patm,Rspes)
    
    T = patm / (rho*Rspes) - 273.15;
    
end

