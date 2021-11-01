%% Path
folders = [genpath('Materials'), ...
    genpath('Numerics'), ...
    genpath('Figures'), ...
    genpath('Scripts'), ...
    genpath('Data')];
addpath(folders);

%% Code
%% Square
square_example;

%% Main Code
project;
draw_stress_and_force_displacement;
draw_plastic_response;