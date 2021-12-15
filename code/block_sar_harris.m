function [key_point_array, key_point_array2, temp_gradient,temp_angle] = block_sar_harris(I_opt,Th_opt_1,Th_opt_2,image_show)
   %% 参数信息
    hist_show=0;
    image_show_point=0;
    sigma = 2;
    ratio = 2^(1/3.);
    Mmax = 1;           % 在单一尺度下的梯度方向和幅值提取
    [M,N] = size(I_opt);
    gradient = zeros(M,N,Mmax);
    angle = zeros(M,N,Mmax);
    HIST_BIN = 18;                                      % 考虑0~180度，分为18个小份。
    SIFT_ORI_PEAK_RATIO=1;                              % 原来为0.9，现在改为1
    key_number=0;             % 关键点数量
    key_number2=0;
    %% 计算单一尺度下的梯度幅值和方向图像，
    i=1;
    scale=sigma*(ratio)^(i-1);
    radius=round(2*scale);
    j=-radius:1:radius;
    k=-radius:1:radius;
    [xarry,yarry]=meshgrid(j,k);

    W = exp(-(abs(xarry)+abs(yarry))/scale);
    W34=zeros(2*radius+1,2*radius+1);
    W12=zeros(2*radius+1,2*radius+1);
    W14=zeros(2*radius+1,2*radius+1);
    W23=zeros(2*radius+1,2*radius+1);

    W34(radius+2:2*radius+1,:)=W(radius+2:2*radius+1,:);
    W12(1:radius,:)=W(1:radius,:);
    W14(:,radius+2:2*radius+1)=W(:,radius+2:2*radius+1);
    W23(:,1:radius)=W(:,1:radius);

    M34=imfilter(I_opt,W34,'replicate');
    M12=imfilter(I_opt,W12,'replicate');
    M14=imfilter(I_opt,W14,'replicate');
    M23=imfilter(I_opt,W23,'replicate');

    Ix=log(M14./M23);
    Iy=log(M34./M12);
    % 得到图像I在两个方向上的梯度
    Ix(find(imag(Ix)))=abs(Ix(find(imag(Ix))));
    Iy(find(imag(Iy)))=abs(Iy(find(imag(Iy))));
    Ix(~isfinite(Ix))=0;
    Iy(~isfinite(Iy))=0;

    temp_gradient=sqrt(Ix.^2+Iy.^2);  
    gradient(:,:,i)=temp_gradient/max(temp_gradient(:));  % 梯度幅值
%     image_3D_output('ROEWA梯度幅值图像', gradient);

    temp_angle=atan(Iy./Ix);
    temp_angle=temp_angle/pi*180;
    temp_angle(temp_angle<0)=temp_angle(temp_angle<0)+180; % 梯度方向
    angle(:,:,i) = temp_angle;
%     image_3D_output('ROEWA梯度方向图像', angle);
    %% 计算相位一致性，根据相位一致性的最大矩和最小矩图进行角点和边缘点的寻找
    [M_pc, m_pc, or,pc,~,~] = phasecong2(I_opt);
    % 对最大矩图和最小矩图进行增强处理
%     a=max(M_pc(:)); b=min(M_pc(:)); M_pc=(M_pc-b)/(a-b);
%     a=max(m_pc(:)); b=min(m_pc(:)); m_pc=(m_pc-b)/(a-b);
    
    temp_gradient=gradient;
    temp_angle=angle;
    opt_harris_function1 = m_pc;    % 最小矩图像，检测PC角点
    opt_harris_function2 = M_pc;    % 最大矩图像，检测PC边缘点
    % 三种图像的立体显现
    if image_show == 1
        image_3D_output('SAR梯度幅值', pc{1});
        image_3D_output('SAR梯度方向', or);
        image_3D_output('SAR的Harris空间', m_pc);
    end

    % 寻找Harris点集,涉及到非极大值抑制
    [M1,N1] = size(opt_harris_function1);
    [M2,N2] = size(opt_harris_function2);
    BORDER_WIDTH = 2;
    length_r = 4;     % 关键点的直径
    mark_r = length_r/2;% 关键点的半径大小
    
    % 返回的关键点的信息
    key_point_array = zeros(M1,6);
    key_point_array2 = zeros(M2,6);
    if image_show_point == 1
        figure,imshow(I_opt),title('原始输入SAR图像');
    end
    temp_current1 = opt_harris_function1;
    temp_current2 = opt_harris_function2;
    
    % 角点
    for i = BORDER_WIDTH:1:M1-BORDER_WIDTH
        for j = BORDER_WIDTH:1:N1-BORDER_WIDTH
            temp = temp_current1(i,j);
            if temp>Th_opt_1 && ... % 与八个邻域值和设定的最大值进行比较，非极大值抑制？
                    temp > temp_current1(i-1,j-1) && temp>temp_current1(i-1,j) && temp>temp_current1(i-1,j+1)...
                    && temp>temp_current1(i,j-1)&& temp>temp_current1(i,j+1)...
                    && temp>temp_current1(i+1,j-1) && temp>temp_current1(i+1,j) && temp>temp_current1(i+1,j+1)
                
                % 根据上述检测出来的特征点，以及梯度幅值，梯度方向图建立特征点描述子，计算各个方向的统计直方图，和主方向的幅值max_value
                [hist,max_value] = calculate_oritation_hist_sar(j,i,2,temp_gradient,temp_angle,HIST_BIN,0);
                % 结合上述坐标建立梯度直方图
                mag_thr=max_value*SIFT_ORI_PEAK_RATIO;  
                for kk=1:1:HIST_BIN
                    % 18个Bin： k1 kk k2
                    if(kk==1)
                        k1=HIST_BIN;
                    else
                        k1=kk-1;
                    end 
                    if(kk==HIST_BIN)
                        k2=1;
                    else
                        k2=kk+1;
                    end
                     if(hist(kk)>hist(k1) && hist(kk)>hist(k2)...
                          && hist(kk)>=mag_thr)
                  
                        bin=kk-1+0.5*(hist(k1)-hist(k2))/(hist(k1)+hist(k2)-2*hist(kk));
                        if(bin<0)
                            bin=HIST_BIN+bin;
                        elseif(bin>=HIST_BIN)
                            bin=bin-HIST_BIN;
                        end
                        key_number=key_number+1;
                        % 将符合条件的点标注出来
                        if image_show_point == 1
                            pos = [j-mark_r,i-mark_r,length_r,length_r];
                            rectangle('Position',pos,'FaceColor',[1,1,0],'curvature',[1,1],'edgecolor','y');
                        end
                        key_point_array(key_number,1)=j;
                        key_point_array(key_number,2)=i;
                        key_point_array(key_number,3)=2;
                        key_point_array(key_number,4)=1;                  % 注意，这里为尺度数。
                        key_point_array(key_number,5)=(180/HIST_BIN)*bin; % 量化方向(0~180度)
                        key_point_array(key_number,6)=hist(kk);           % 主方向幅值
                    end
                end
            end
        end
    end
    fprintf("SAR图像角点数 %d 个。\n",key_number);
    
    % 边缘点
    for i = BORDER_WIDTH:1:M2-BORDER_WIDTH
        for j = BORDER_WIDTH:1:N2-BORDER_WIDTH
            temp = temp_current2(i,j);
            if temp>Th_opt_2 && ... % 与八个邻域值和设定的最大值进行比较，非极大值抑制？
                    temp > temp_current2(i-1,j-1) && temp>temp_current2(i-1,j) && temp>temp_current2(i-1,j+1)...
                    && temp>temp_current2(i,j-1)&& temp>temp_current2(i,j+1)...
                    && temp>temp_current2(i+1,j-1) && temp>temp_current2(i+1,j) && temp>temp_current2(i+1,j+1)
                
                % 根据上述检测出来的特征点，以及梯度幅值，梯度方向图建立特征点描述子，计算各个方向的统计直方图，和主方向的幅值max_value
                [hist,max_value] = calculate_oritation_hist_sar(j,i,2,temp_gradient,temp_angle,HIST_BIN,0);
                % 结合上述坐标建立梯度直方图
                mag_thr=max_value*SIFT_ORI_PEAK_RATIO;  
                for kk=1:1:HIST_BIN
                    % 18个Bin： k1 kk k2
                    if(kk==1)
                        k1=HIST_BIN;
                    else
                        k1=kk-1;
                    end 
                    if(kk==HIST_BIN)
                        k2=1;
                    else
                        k2=kk+1;
                    end
                     if(hist(kk)>hist(k1) && hist(kk)>hist(k2)...
                          && hist(kk)>=mag_thr)
                  
                        bin=kk-1+0.5*(hist(k1)-hist(k2))/(hist(k1)+hist(k2)-2*hist(kk));
                        if(bin<0)
                            bin=HIST_BIN+bin;
                        elseif(bin>=HIST_BIN)
                            bin=bin-HIST_BIN;
                        end
                        key_number2=key_number2+1;
                        % 将符合条件的点标注出来
                        if image_show_point == 1
                            pos = [j-mark_r,i-mark_r,length_r,length_r];
                            rectangle('Position',pos,'FaceColor',[1,0,0],'curvature',[1,1],'edgecolor','r');
                        end
                        key_point_array2(key_number2,1)=j;
                        key_point_array2(key_number2,2)=i;
                        key_point_array2(key_number2,3)=2;
                        key_point_array2(key_number2,4)=1;                  % 注意，这里为尺度数。
                        key_point_array2(key_number2,5)=(180/HIST_BIN)*bin; % 量化方向(0~180度)
                        key_point_array2(key_number2,6)=hist(kk);           % 主方向幅值
                    end
                end
            end
        end
    end
    fprintf("SAR图像边缘点数 %d 个。\n",key_number2);
    key_point_array2=key_point_array2(1:key_number2,1:6);
    
    fprintf("SAR图像特征点数 %d 个。\n",key_number + key_number2);
    % 得到key_point_array后如何去掉重复点和相应的重复直方图特征
%     key_point_array=key_point_array(1:key_number,1:6);
%     key_point_array=unique(key_point_array, 'rows');
%     fprintf("SAR图像去掉重复点后特征点数为 %d 个。\n",size(key_point_array, 1));
end

