%% Main
clear; close all; clc

%% Parameters

f = 0.105; % camera focal length (m)
p = 5.95e-6; % pixel pitch (m)
ZB = 1.09; % distance between background and focal point (m) 
Lx = 0.400; % paper size x (m)
Ly = 0.400; % paper size y (m)
dot_size_pixels = 5; % dot size (px)
dot_spacing = 3; % dot spacing (dots)

%% Processing

generate_background(f,p,ZB,Lx,Ly,dot_size_pixels,dot_spacing)