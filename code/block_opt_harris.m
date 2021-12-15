function [key_point_array, key_point_array2,temp_gradient,temp_angle] = block_opt_harris(I_opt,Th_opt_1,Th_opt_2,image_show)
% 同时使用M_pc和m_pc
% 阈值Th_opt应该改一改，使得，Th_opt_1,Th_opt_2
   %% 参数信息
    hist_show=0;
    image_show_point=0; % 特征点检测结果显示设置：0不显示，1显示
    sigma = 2;
    ratio = 2^(1/3.);
    [M,N] = size(I_opt);
    Mmax = 1;
    gradient = zeros(M,N,Mmax);
    angle = zeros(M,N,Mmax);
    HIST_BIN = 18;                                      % 考虑0~180度，分为18个小份。
    SIFT_ORI_PEAK_RATIO=1;
    key_number=0;             % 关键点数量
    key_number2=0;
    %% 计算单一尺度下的梯度幅值和方向图
    ii = 1;
    scale = sigma*(ratio)^(ii-1);                                        % 后期需要调整
    radius=round(scale);
    jj=-radius:1:radius;
    k=-radius:1:radius;
    [xarry,yarry]=meshgrid(jj,k);
    W=exp(-((xarry.*xarry)+(yarry.*yarry))/(2*scale));
    W2=zeros(2*radius+1,2*radius+1);
    W1=zeros(2*radius+1,2*radius+1);
    W2(radius+2:2*radius+1,:)=W(radius+2:2*radius+1,:);
    W2(1:radius,:)=-W(1:radius,:);
    W1(:,radius+2:2*radius+1)=W(:,radius+2:2*radius+1);
    W1(:,1:radius)=-W(:,1:radius);

    Ix=imfilter(I_opt,W1,'replicate');
    Iy=imfilter(I_opt,W2,'replicate');

    Ix(find(imag(Ix)))=abs(Ix(find(imag(Ix))));
    Iy(find(imag(Iy)))=abs(Iy(find(imag(Iy))));
    Ix(~isfinite(Ix))=0;
    Iy(~isfinite(Iy))=0;

    temp_gradient=sqrt(Ix.^2+Iy.^2);      
    gradient(:,:,ii)=temp_gradient/max(temp_gradient(:));      % 梯度幅值图像
%     image_3D_output('Sobel梯度幅值图像', gradient);

    temp_angle=atan(Iy./Ix);
    temp_angle=temp_angle/pi*180;
    temp_angle(temp_angle<0)=temp_angle(temp_angle<0)+180;  % 梯度方向图像
    angle(:,:,ii)=temp_angle;
%     image_3D_output('Sobel梯度方向图像', angle);
    %% 在PC最大矩和最小矩图上寻找极值点
    [M_pc, m_pc, or,pc,~,~] = phasecong2(I_opt);
    % 对最大矩图和最小矩图进行增强处理
%     a=max(M_pc(:)); b=min(M_pc(:)); M_pc=(M_pc-b)/(a-b);
%     a=max(m_pc(:)); b=min(m_pc(:)); m_pc=(m_pc-b)/(a-b);
    
    temp_gradient=gradient;
    temp_angle = angle;
    opt_harris_function1 = m_pc;   % 最小矩图 
    opt_harris_function2 = M_pc;   % 最大矩图
    % 三种图像的立体显现                  
    if image_show == 1
        image_3D_output('光学梯度幅值', pc);
        image_3D_output('光学梯度方向', or);
        image_3D_output('光学Harris空间', m_pc);
    end

    % 寻找Harris点集,涉及到非极大值抑制
    [M1,N1] = size(opt_harris_function1);
    [M2,N2] = size(opt_harris_function2);
    BORDER_WIDTH = 2;
    length_r = 4;     % 关键点的直径
    mark_r = length_r/2;% 关键点的半径大小
    
    % 返回的关键点的信息
    key_point_array = zeros(M1,6);        % 
    key_point_array2 = zeros(M2,6);

    if image_show_point == 1
        figure,imshow(I_opt),title('原始输入光学图像');
    end
    temp_current1 = opt_harris_function1;         % 取得最小矩图
    temp_current2 = opt_harris_function2;         % 取得最大矩图
    
    % 在最小矩图上寻找极值点，角点颜色为黄色
    for ii = BORDER_WIDTH:1:M1-BORDER_WIDTH
        for jj = BORDER_WIDTH:1:N1-BORDER_WIDTH
            temp = temp_current1(ii,jj);
            if temp>Th_opt_1 && ... % 与八个邻域值和设定的最大值进行比较，非极大值抑制？
                    temp > temp_current1(ii-1,jj-1) && temp>temp_current1(ii-1,jj) && temp>temp_current1(ii-1,jj+1)...
                    && temp>temp_current1(ii,jj-1)&& temp>temp_current1(ii,jj+1)...
                    && temp>temp_current1(ii+1,jj-1) && temp>temp_current1(ii+1,jj) && temp>temp_current1(ii+1,jj+1)
                
                % 根据上述检测出来的特征点，以及梯度幅值，梯度方向图建立特征点描述子，计算各个方向的统计直方图，和主方向的幅值max_value
                if ii == BORDER_WIDTH+10 && jj == BORDER_WIDTH+10
                    hist_show = 1;
                else
                    hist_show  = 0;
                end
                [hist,max_value] = calculate_oritation_hist_sar(jj,ii,2,temp_gradient,temp_angle,HIST_BIN,0);
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
                            pos = [jj-mark_r,ii-mark_r,length_r,length_r];
                            rectangle('Position',pos,'FaceColor',[1,1,0],'curvature',[1,1],'edgecolor','y');
                        end
                        key_point_array(key_number,1)=jj;
                        key_point_array(key_number,2)=ii;
                        key_point_array(key_number,3)=2;
                        key_point_array(key_number,4)=1;                  % 注意，这里为尺度数。
                        key_point_array(key_number,5)=(180/HIST_BIN)*bin; % 量化方向(0~180度)
                        key_point_array(key_number,6)=hist(kk);           % 主方向幅值
                    end
                end
            end
        end
    end
    fprintf("光学图像角点数 %d 个。\n",key_number);
    key_point_array=key_point_array(1:key_number,1:6);
    
    % 在最大矩图上寻找极值点，边缘点颜色为红色
    for ii = BORDER_WIDTH:1:M2-BORDER_WIDTH
        for jj = BORDER_WIDTH:1:N2-BORDER_WIDTH
            temp = temp_current2(ii,jj);
            if temp>Th_opt_2 && ... % 与八个邻域值和设定的最大值进行比较，非极大值抑制？
                    temp > temp_current2(ii-1,jj-1) && temp>temp_current2(ii-1,jj) && temp>temp_current2(ii-1,jj+1)...
                    && temp>temp_current2(ii,jj-1)&& temp>temp_current2(ii,jj+1)...
                    && temp>temp_current2(ii+1,jj-1) && temp>temp_current2(ii+1,jj) && temp>temp_current2(ii+1,jj+1)
                
                % 根据上述检测出来的特征点，以及梯度幅值，梯度方向图建立特征点描述子，计算各个方向的统计直方图，和主方向的幅值max_value
                [hist,max_value] = calculate_oritation_hist_sar(jj,ii,2,temp_gradient,temp_angle,HIST_BIN,0);
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
                            pos = [jj-mark_r,ii-mark_r,length_r,length_r];
                            rectangle('Position',pos,'FaceColor',[1,0,0],'curvature',[1,1],'edgecolor','r');
                        end
                        key_point_array2(key_number2,1)=jj;
                        key_point_array2(key_number2,2)=ii;
                        key_point_array2(key_number2,3)=2;
                        key_point_array2(key_number2,4)=1;                  % 注意，这里为尺度数。
                        key_point_array2(key_number2,5)=(180/HIST_BIN)*bin; % 量化方向(0~180度)
                        key_point_array2(key_number2,6)=hist(kk);           % 主方向幅值
                    end
                end
            end
        end
    end
    fprintf("光学图像边缘点数 %d 个。\n",key_number2);
    key_point_array2=key_point_array2(1:key_number2,1:6);
    
    fprintf("光学图像特征点数 %d 个。\n",key_number+key_number2);
    % 得到key_point_array去掉重复点和相应的重复直方图特征
%     key_point_array=key_point_array(1:key_number,1:6);
%     key_point_array=unique(key_point_array, 'rows');
%     fprintf("光学图像去掉重复点后特征点数为 %d 个。\n",size(key_point_array, 1));
end

