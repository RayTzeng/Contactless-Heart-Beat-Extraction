%Demo

clc,clear; close all
addpath(genpath([cd '\tools\'])) %additional required function

DataDirectory         = [cd '\testvideo\']; %default determined off this script's directory
VideoFile             = [DataDirectory 'run_10min_3.MOV'];% Video file path.
FS                    = 30; % Video Framerate(fps)
StartTime             = 0;% Timepoint at which to start process (default = 0 seconds).
Duration              = 600;% Duration of the time window to process (default = 60 seconds).
PlotTF                = true;% Boolean to turn plotting results on or off.
%Histogram             = ['img']

[BVP,BVP_I,PR,M,Sx,W] = CrPPG(VideoFile, FS, StartTime, Duration, PlotTF); 
