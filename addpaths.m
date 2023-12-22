clear all;
close all;
clc;

%%-----------------------------------------
% Add libraries
%------------------------------------------

% Tensorlab lib.
addpath("libraries/tensor_toolbox-v3.4/");

% Tensorbox lib.
addpath("libraries/TENSORBOX/");

addpath("libraries/NN_Comp_Toolbox/NN_Comp_Toolbox/")

addpath("libraries/NN_Comp_Toolbox_2/NN_Comp_Toolbox/")

addpath("libraries/NNtensorPackage_mfiles/NNtensorPackage_mfiles/")

% Util functions
addpath("Utils/");
addpath("Datasets/");
% addpath(genpath(pwd))

% Data set
addpath("Datasets/AVIRIS_Cuprite/")