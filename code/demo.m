clear all;close all;clc;
use_gamma = 0;        % 1：使用gamma增强；0：不使用gamma增强
%% 正常匹配图像对参数
I_opt = imread('..\data\P3_1.tif');
I_sar = imread('..\data\P3_2.tif'); 
Th_opt1 = 0.0001;             
Th_opt2 = 0.0001;              
Th_sar1 = 0.0001;             
Th_sar2 = 0.0001; 

%% 不同旋转角度差异图像对匹配
% I_opt = imread('..\data\P3_1.tif');
% I_opt=imread('..\data\P3_1_18度旋转.tif');
% I_opt=imread('..\data\P3_1_36度旋转.tif');
% I_opt=imread('..\data\P3_1_54度旋转.tif');
% I_opt=imread('..\data\P3_1_72度旋转.tif');
% I_sar = imread('..\data\P3_2_small.tif'); 
% 除了36度都适用
% Th_opt1 = 0.0001;               
% Th_opt2 = 0.0001;               
% Th_sar1 = 0.0001;             
% Th_sar2 = 0.0001; 
% 36度旋转图像参数
% Th_opt1 = 0.01;               
% Th_opt2 = 0.01;              
% Th_sar1 = 0.00001;            
% Th_sar2 = 0.00001; 

%% 不同尺度比图像对匹配匹配
% I_opt = imread('..\data_rot\P3_1.tif');      % 1.0
% I_opt=imread('..\data_scale\P3_1_450.tif');  % 0.9 
% I_opt=imread('..\data_scale\P3_1_400.tif');  % 0.8
% I_opt=imread('..\data_scale\P3_1_350.tif');  % 0.7
% I_sar = imread('..\data_scale\P3_2.tif');    % 1.0
% I_sar = imread('..\data_scale\P3_2_450.tif');  % 1.1
% I_sar = imread('..\data_scale\P3_2_400.tif');  % 1.2
% 500~350
% Th_opt1 = 0.0001;               % 光学阈值0.05，0.01
% Th_opt2 = 0.0001;               % 调好阈值进行
% Th_sar1 = 0.0001;               % SAR阈值0.01，0.005
% Th_sar2 = 0.0001; 
% 300
% Th_opt1 = 0.000001;               % 光学阈值0.05，0.01
% Th_opt2 = 0.000001;               % 调好阈值进行
% Th_sar1 = 0.000001;               % SAR阈值0.01，0.005
% Th_sar2 = 0.000001; 


%% 初始图像预处理
image_show = 0;
is_keypoints_refine = false;
if size(I_opt,3)>1
    I_opt = rgb2gray(I_opt);   
end
if size(I_sar,3)>1
    I_sar = rgb2gray(I_sar);
end
I_opt = im2double(I_opt);
I_sar = im2double(I_sar);

I_opt = I_opt + 0.001;       
I_sar = I_sar + 0.001; 

t1 = clock;
%% PC最大矩点和最小矩点检测
[key_point_opt, key_point_opt2,opt_gradient,opt_angle] = block_opt_harris(I_opt,Th_opt1,Th_opt2,image_show);
[key_point_sar, key_point_sar2,sar_gradient,sar_angle] = block_sar_harris(I_sar,Th_sar1,Th_sar2,image_show);
%% 特征描述算法，使用ROWEWA-GLOH和Sobel-GLOH分别对SAR图像和光学图像的特征点建立描述符
% 对最小矩点建立特征描述符
[describe_1, locs_1] = build_describe(key_point_opt,opt_gradient,opt_angle);        % 光学图像
[describe_2, locs_2] = build_describe(key_point_sar,sar_gradient,sar_angle);        % SAR图像
% 对最大矩点建立特征描述符
[describe_1_2, locs_1_2] = build_describe(key_point_opt2,opt_gradient,opt_angle);   % 光学图像
[describe_2_2, locs_2_2] = build_describe(key_point_sar2,sar_gradient,sar_angle);   % SAR图像
%% 特征匹配+尝试加入GMS对NNDR后的初始匹配点对进行筛选
[solution,cor1,cor2,cor1_2,cor2_2,cor11,cor22,cor11_2,cor22_2] = match(I_opt,I_sar,describe_1,locs_1,describe_2,locs_2,...
                                                     describe_1_2, locs_1_2, describe_2_2, locs_2_2, false);    
t2 = clock;
fprintf('总运行时间为 %f.\n', etime(t2,t1));
%% 图像拼接与融合
if use_gamma==1            % 因SAR图像总体较暗，在最终结果中可以使用伽马校正，提高SAR图像亮度
    I_sar = sqrt(I_sar);
end
image_fusion(I_sar,I_opt,solution);
