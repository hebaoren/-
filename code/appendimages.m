function  [im3]=appendimages(image1, image2,correspond1,correspond2,correspond1_2,correspond2_2)

rows1 = size(image1,1);
rows2 = size(image2,1);

col1=size(image1,2);
col2=size(image2,2);

if (rows1 < rows2)
     image1(rows1+1:rows2,1:col1,:) = 0;
elseif(rows1 >rows2)
     image2(rows2+1:rows1,1:col2,:) = 0;
end

temp1=size(image1,3);
temp2=size(image2,3);
if(temp1==1 && temp2==3)
    image2=rgb2gray(image2);
elseif(temp1==3 && temp2==1)
    image1=rgb2gray(image1);
end
im3 = [image1 image2]; 

% imwrite(im3,'..\data\SAR_OPT_fusion.png');   % 存储融合后的图片

fhand=figure;
imshow(im3,'border','tight','initialmagnification','fit');
title(['left is the reference --- the number of pairs ',num2str(size(correspond1,1)),' --- right is the To be registered']);
set (gcf,'Position',[0,0,size(im3,2) size(im3,1)]);
axis normal;

hold on;
cols1 = size(image1,2);

num_point=1;
num_line=1;
num_point2=2;
num_line2=2;

% 画线，注意这里的蓝线
for i = 1: size(correspond1,1)
    if(num_line==1)%red
        line([correspond1(i,1) correspond2(i,1)+cols1], ...
             [correspond1(i,2) correspond2(i,2)], 'Color', 'r','LineWidth',1);
    elseif(num_line==2)%green
        line([correspond1(i,1) correspond2(i,1)+cols1], ...
             [correspond1(i,2) correspond2(i,2)], 'Color', 'g','LineWidth',1); 
    elseif(num_line==3)%blue
        line([correspond1(i,1) correspond2(i,1)+cols1], ...
             [correspond1(i,2) correspond2(i,2)], 'Color', 'b','LineWidth',1); 
    end
end
for i = 1: size(correspond1_2,1)
    if(num_line2==1)%red
        line([correspond1_2(i,1) correspond2_2(i,1)+cols1], ...
             [correspond1_2(i,2) correspond2_2(i,2)], 'Color', 'r','LineWidth',1);
    elseif(num_line2==2)%green
        line([correspond1_2(i,1) correspond2_2(i,1)+cols1], ...
             [correspond1_2(i,2) correspond2_2(i,2)], 'Color', 'y','LineWidth',1); 
    elseif(num_line2==3)%blue
        line([correspond1_2(i,1) correspond2_2(i,1)+cols1], ...
             [correspond1_2(i,2) correspond2_2(i,2)], 'Color', 'b','LineWidth',1); 
    end
end

% 标注关键点
% 首先标注最小矩图的角点
for i = 1: size(correspond1,1)
    length_r = 2;     % 关键点的直径
    mark_r = length_r/2;% 关键点的半径大小
    pos1 = [correspond1(i,1)-mark_r,correspond1(i,2)-mark_r,length_r,length_r];
    pos2 = [correspond2(i,1)-mark_r+cols1, correspond2(i,2)-mark_r,length_r,length_r];
    if num_point == 1
        rectangle('Position',pos1,'FaceColor',[1,0,0],'curvature',[1,1],'edgecolor','r');
        rectangle('Position',pos2,'FaceColor',[1,0,0],'curvature',[1,1],'edgecolor','r');
    end
    if num_point == 3
        rectangle('Position',pos1,'curvature',[1,1],'edgecolor','g');
        rectangle('Position',pos2,'curvature',[1,1],'edgecolor','g');
    end
    if num_point == 2
        rectangle('Position',pos1,'curvature',[1,1],'edgecolor','y');
        rectangle('Position',pos2,'curvature',[1,1],'edgecolor','y');
    end
end
% 然后标注最大矩图的边缘点特征，通过相应的手段进行点对滤除
for i = 1: size(correspond1_2,1)
    length_r = 2;     % 关键点的直径
    mark_r = length_r/2;% 关键点的半径大小
    pos1 = [correspond1_2(i,1)-mark_r,correspond1_2(i,2)-mark_r,length_r,length_r];
    pos2 = [correspond2_2(i,1)-mark_r+cols1, correspond2_2(i,2)-mark_r,length_r,length_r];
    if num_point2 == 1
        rectangle('Position',pos1,'curvature',[1,1],'edgecolor','y');
        rectangle('Position',pos2,'curvature',[1,1],'edgecolor','y');
    end
    if num_point2 == 2
        rectangle('Position',pos1,'FaceColor',[1,1,0],'curvature',[1,1],'edgecolor','y');
        rectangle('Position',pos2,'FaceColor',[1,1,0],'curvature',[1,1],'edgecolor','y');
    end
    if num_point2 == 2
        rectangle('Position',pos1,'curvature',[1,1],'edgecolor','b');
        rectangle('Position',pos2,'curvature',[1,1],'edgecolor','b');
    end
end
hold off;
end






