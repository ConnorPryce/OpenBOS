%% Main
clear; close all; clc

%% Parameters

f = 0.105; % camera focal length (m)
P = 5.95e-6; % pixel pitch (m)
ZA = 0.53; % distance between object and focal point (m)
ZD = 0.56; % distance between object and background (m)
D = 0.005; % medium depth (m)
Lref = 0.057; % reference length for calibration (m)
n0 = 1.000273; % air index of refraction (-)
G = 0.2228e-3; % Gladstone-dale constant for air (m^3/kg)
patm = 101e3; % atmospheric air pressure (Pa)
R = 287; % specific gas constant for air (J/kgK)
window_size = 31; % OF window size (px)
overlap = 25; % OF window overlap (px)

%% Load Video

[frames,~] = load_video("candle.MOV");
image_mask = create_mask(frames);

I = frames(:,:,2:end);
Iref = frames(:,:,1);
IrefL = frames(:,:,2);

M = calculate_magnification_factor(IrefL,Lref);

%% Processing

[X,Y,ds,gx,gy,domain_mask] = pixel_displacements(Iref,I,window_size,overlap,image_mask);
[cx,cy,~] = image_registration(X,Y,gx,gy,domain_mask);
[dndx,dndy] = refraction_gradients(gx+cx,gy+cy,n0,f,ZA,ZD,D,P);
[Xt,Yt,dst] = transform_coordinates(M,X,Y,ds);
n = poisson_reconstruction(dst,dndx,dndy,domain_mask);
n = shift_solution(n,n0,8,1);
rho = calculate_density(n,G);
T = calculate_temperature(rho,patm,R);

%% Plotting 

xmin = min(Xt(:)); xmax = max(Xt(:));
ymin = min(Yt(:)); ymax = max(Yt(:));
rhomin = 0.8; rhomax = 1.3;
Tmin = 15; Tmax = 200;

figure(1)

for i = 1:size(n,3)

tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
set(gcf,'Color','white')
set(gcf, 'Position',  [100, 100, 500*2, 400*2])

nexttile()
imshow(frames(:,:,i+1))
title("Input Video",Interpreter="latex")

nexttile()
surf(Xt,Yt,hypot(gx(:,:,i)+cx(:,:,i),gy(:,:,i)+cy(:,:,i))); hold on
ax = gca;
ax.XColor = 'none';
ax.YColor = 'none';
axis equal 
view(2)
shading interp
grid off
xlim([xmin,xmax])
ylim([ymin,ymax])
colormap(gca,parula)
title("Displacement Magnitude $$(m)$$",Interpreter="latex")

nexttile()
surf(Xt,Yt,rho(:,:,i)); hold on
ax = gca;
ax.XColor = 'none';
ax.YColor = 'none';
axis equal 
view(2)
shading interp
grid off
xlim([xmin,xmax])
ylim([ymin,ymax])
clim([rhomin,rhomax])
colorbar
colormap(gca,parula)
title("Density $$(kg/m^3)$$",Interpreter="latex")

nexttile()
surf(Xt,Yt,T(:,:,i)); hold on
ax = gca;
ax.XColor = 'none';
ax.YColor = 'none';
axis equal 
view(2)
shading interp
grid off
xlim([xmin,xmax])
ylim([ymin,ymax])
clim([Tmin,Tmax])
colorbar
colormap(gca,hot)
title("Temperature $$(^\circ C)$$",Interpreter="latex")

pause(0.0001)

end 


