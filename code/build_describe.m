function [descriptors, locs] = build_describe(key_point_sar,sar_gradient,sar_angle)
% sar_describe：SAR图像中的特征点的描述子建立
%   输入:特征点坐标,光学梯度图像,光学角度图像
    circle_bin=8;
    LOG_DESC_HIST_BINS=8;
    M = size(key_point_sar, 1);
    d = circle_bin;
    n = LOG_DESC_HIST_BINS;
    % 为描述子分配空间
    descriptors = zeros(M,(2*d+1)*n);
    locs = key_point_sar;
    
    % 计算各个特征点的描述子，并将其存储到一个矩阵下进行输出
    parfor i=1:M
        x = key_point_sar(i,1);
        y = key_point_sar(i,2);
        scale = key_point_sar(i,3);
        layer = key_point_sar(i,4);
        main_angle = key_point_sar(i,5);
        current_gradient = sar_gradient;
        % 注意max(a(:))和max(max(a))的结果一样，max(a)输出为矩阵的某一行
        current_gradient=current_gradient/max(current_gradient(:));  % opt_angle(:,:,layer)
        current_angle=sar_angle;                                  % sar_angle(:,:,layer)
        
%         descriptors(i,1:(2*d+1)*n)=calc_log_polar_descriptor(current_gradient,current_angle,...
%             x,y,scale,main_angle,d,n);
        descriptors(i,:)=calc_log_polar_descriptor(current_gradient,current_angle,...
            x,y,scale,main_angle,d,n);
    end
end

