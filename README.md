# TRT_Filter
Thresholding and Regularization Techniques for Image Denoising

GUI requires MATLAB 2021b

===========================================================


% For the first time you run the app, you need to set up MEX enviroment 

% by running the following command:

% install

% A mex file will be generated, usually tv_restore.mexmaci64,

% tv_restore.mexwindow64, etc.


% =============================

% If you want to test the algorithm without GUI, run the below code:


% clc;

% close all;

% I = imread('232038.jpeg');

% [In, Id] = TRTDenoise(I, .7);

% imshow([I, In, Id]);


% If you want to use GUI, run the below code:

TRTapp
