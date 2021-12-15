%% 增加纹理欠丰富的SAR和光学图像对匹配实验
% 
clear all;close all;clc;
use_gamma = 0;        % 1：使用gamma增强；0：不使用gamma增强
%% 不同的光学和SAR图像对
% I_opt = imread('..\data\1山区-光学.jpg');
% I_sar = imread('..\data\1山区-SAR.jpg');

% I_opt = imread('..\data\OPT_Farmland_PSO.jpg');
% I_sar = imread('..\data\SAR_Farmland_PSO.jpg');

% I_sar = imread('..\data\SAR_sen.png');
% I_opt = imread('..\data\optical_ref.png');

% I_opt = imread('..\data\t2.tif');
% I_sar = imread('..\data\t1.tif');
% use_gamma = 1;

% I_opt = imread('..\data\P3_1_small.tif');
% I_sar = imread('..\data\P3_2_small.tif'); 

% I_opt = imread('..\data\OPT_2.png');
% I_sar = imread('..\data\SAR_2.png'); 
% I_opt = imread('..\data\光学1.png');
% I_sar = imread('..\data\sar1.png'); 

% 欠纹理图像对 P-E图像对
I_opt = imread('..\data\OPT_2 - 副本.png');
I_sar = imread('..\data\SAR_2 - 副本.png');  
% shi'yong'ruo'wen'li'qu
%% 相关参数设置和图像预处理
% 针对t1和t2图像对的最佳参数，75个，RMSE=1.653485 and 1.692323，效果比OS-SIFT要好
% Th_opt1 = 0.008;               % 光学阈值0.05，0.01
% Th_opt2 = 0.0001;               % 调好阈值进行
% Th_sar1 = 0.004;               % SAR阈值0.01，0.005
% Th_sar2 = 0.001; 

% 针对t1_small和t2_small
% Th_opt1 = 0.0002;               % 光学阈值0.05，0.01
% Th_opt2 = 0.00000001;               % 调好阈值进行
% Th_sar1 = 0.0005;               % SAR阈值0.01，0.005
% Th_sar2 = 0.0005; 

% 针对山区图像对的最佳参数
% Th_opt1 = 0.003;  
% Th_opt2 = 0.002;
% Th_sar1 = 0.0001;
% Th_sar2 = 0.00001;

% 针对农田图像对的最佳参数
% Th_opt1 = 0.001;
% Th_opt2 = 0.001;
% Th_sar1 = 0.0001;
% Th_sar2 = 0.001;

% % 针对P3_1_small和P3_2_small图像的最佳参数，对原图像P3_1和P_3_2的匹配结果也很好
% Th_opt1 = 0.0001;               % 光学阈值0.05，0.01
% Th_opt2 = 0.0001;               % 调好阈值进行
% Th_sar1 = 0.0001;               % SAR阈值0.01，0.005
% Th_sar2 = 0.0001; 

% 补充弱纹理实验
% 对OPT_2和SAR_2图像对的最佳参数，效果比OS-SIFT要好
% Th_opt1 = 0.002;               % 光学阈值0.05，0.01
% Th_opt2 = 0.001;               % 调好阈值进行
% Th_sar1 = 0.002;               % SAR阈值0.01，0.005
% Th_sar2 = 0.001; 

% 针对农P-E图像对对的最佳参数
Th_opt1 = 0.002;               % 光学阈值0.05，0.01
Th_opt2 = 0.001;               % 调好阈值进行
Th_sar1 = 0.002;               % SAR阈值0.01，0.005
Th_sar2 = 0.001; 

image_show = 0;
is_keypoints_refine = false;  % 试着精调一下
% 预处理
if size(I_opt,3)>1
    I_opt = rgb2gray(I_opt);   
end
if size(I_sar,3)>1
    I_sar = rgb2gray(I_sar);
end
I_opt = im2double(I_opt);
I_sar = im2double(I_sar);

I_opt = I_opt + 0.001;       
I_sar = I_sar + 0.001;       % 防止相应值为0，弱纹理区域

t1 = clock;
%% PC角点和边缘点检测算法
[key_point_opt, key_point_opt2,opt_gradient,opt_angle] = block_opt_harris(I_opt,Th_opt1,Th_opt2,image_show);
[key_point_sar, key_point_sar2,sar_gradient,sar_angle] = block_sar_harris(I_sar,Th_sar1,Th_sar2,image_show);
%% 特征描述算法，使用ROWEWA-GLOH和Sobel-GLOH分别对SAR图像和光学图像的特征建立描述符
[describe_1, locs_1] = build_describe(key_point_opt,opt_gradient,opt_angle);        % 光学图像
[describe_2, locs_2] = build_describe(key_point_sar,sar_gradient,sar_angle);        % SAR图像
% 建立边缘点描述子
[describe_1_2, locs_1_2] = build_describe(key_point_opt2,opt_gradient,opt_angle);   % 光学图像
[describe_2_2, locs_2_2] = build_describe(key_point_sar2,sar_gradient,sar_angle);   % SAR图像
%% 特征匹配+尝试加入GMS对NNDR后的初始匹配点对进行筛选
[solution,cor1,cor2,cor1_2,cor2_2,cor11,cor22,cor11_2,cor22_2] = match(I_opt,I_sar,describe_1,locs_1,describe_2,locs_2,...
                                                     describe_1_2, locs_1_2, describe_2_2, locs_2_2, false);    
t2 = clock;
fprintf('总运行时间为 %f.\n', etime(t2,t1));
%% 图像拼接与融合
if use_gamma==1
    I_sar = sqrt(I_sar);
end
image_fusion(I_sar,I_opt,solution);
