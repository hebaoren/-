function [solution_sum,cor1,cor2,cor1_2,cor2_2,cor11,cor22,cor11_2,cor22_2]= match(im1, im2,des1,loc1,des2,loc2, ...
                                                          des1_2,loc1_2,des2_2,loc2_2,is_multi_region)
distRatio = 0.9;   % 可以适当增加阈值，已得到相应的改进结果
if is_multi_region==false
    des2t = des2';
    for i = 1 : size(des1,1)
        dotprods = des1(i,:) * des2t;
        [vals,indx] = sort(acos(dotprods));
        if (vals(1) < distRatio * vals(2))
            match(i) = indx(1);
        else
            match(i) = 0;
        end
    end
    num = sum(match > 0);
    [~,point1,point2]=find(match);
    cor11=loc1(point1,[1 2 3 4 5]);cor22=loc2(point2,[1 2 3 4 5]);         % droupt values
    cor11=[cor11 point2'];cor22=[cor22 point2'];
else
    des1x=des1(:,1:136);des1y=des1(:,1+136:136+136);des1z=des1(:,1+136*2:136+136*2);
    des2x=des2(:,1:136);des2y=des2(:,1+136:136+136);des2z=des2(:,1+136*2:136+136*2);
    for i = 1 : size(des1,1)
        dx = des1x(i,:) * des2x';dy = des1y(i,:) * des2y';dz = des1z(i,:) * des2z';
        dotprods = (dx+dz+dy)/3;
        [vals,indx] = sort(acos(dotprods));%vals为升序排序的元素，indx为原来的位置
        if (vals(1) < distRatio * vals(2))
            match(i) = indx(1);
        else
            match(i) = 0;
        end
    end
%     disp(size(match));
    num = sum(match > 0);
    [~,point1,point2]=find(match);
    cor11=loc1(point1,[1 2 3 4 5]);cor22=loc2(point2,[1 2 3 4 5]);
    cor11=[cor11 point2'];cor22=[cor22 point2'];
end

%% 
if is_multi_region == false
    des2t = des2_2';
    for j = 1 : size(des1_2,1)
        dotprods = des1_2(j,:) * des2t;
        [vals,indx] = sort(acos(dotprods));
        if (vals(1) < distRatio * vals(2))
            match(j) = indx(1);
        else
            match(j) = 0;
        end
    end
    num = sum(match > 0);
    [~,point1,point2]=find(match);
    cor11_2=loc1_2(point1,[1 2 3 4 5]);cor22_2=loc2_2(point2,[1 2 3 4 5]);
    cor11_2=[cor11_2 point2'];cor22_2=[cor22_2 point2'];
else
    des1x=des1_2(:,1:136);des1y=des1_2(:,1+136:136+136);des1z=des1_2(:,1+136*2:136+136*2);
    des2x=des2(:,1:136);des2y=des2(:,1+136:136+136);des2z=des2(:,1+136*2:136+136*2);
    for j = 1 : size(des1_2,1)
        dx = des1x(j,:) * des2x';dy = des1y(j,:) * des2y';dz = des1z(j,:) * des2z';
        dotprods = (dx+dz+dy)/3;
        [vals,indx] = sort(acos(dotprods));%vals为升序排序的元素，indx为原来的位置
%         if i==1
%             disp(size(dotprods));
%             disp(size(vals));
%             disp(size(indx));
%         end
        if (vals(1) < distRatio * vals(2))
            match(j) = indx(1);
        else
            match(j) = 0;
        end
    end
%     disp(size(match));
    num = sum(match > 0);
    [~,point1,point2]=find(match);
    cor11_2=loc1(point1,[1 2 3 4 5]);cor22_2=loc2(point2,[1 2 3 4 5]);
    cor11_2=[cor11_2 point2'];cor22_2=[cor22_2 point2'];
end

%% 联合匹配点
% 注意，这里输入的是两种匹配点，前者为角点，后者为边缘特征
appendimages(im2,im1,cor22,cor11,cor22_2,cor11_2);

% Remove duplicate points
uni1=[cor11(:,[1 2]),cor22(:,[1 2])];
[~,i,~]=unique(uni1,'rows','first');
cor11=cor11(sort(i)',:);
cor22=cor22(sort(i)',:);

uni1=[cor11_2(:,[1 2]),cor22_2(:,[1 2])];
[~,j,~]=unique(uni1,'rows','first');
cor11_2=cor11_2(sort(j)',:);
cor22_2=cor22_2(sort(j)',:);

num1 = size(cor11_2,1) + size(cor11,1);
fprintf('NNDR found %d matchs.\n', num1);
button=appendimages(im2,im1,cor22,cor11,cor22_2,cor11_2);

%% GMS
% obj = py.importlib.import_module('gms_matcher');  % 
% py.importlib.reload(obj);
% py.gms_matcher_save.orb_gms_matcher();

%% RANSAC
% [solution, inliers] = ransac_point_fit(cor11(:,1:2)', cor22(:,1:2)', 0.01);rmse=0;
% [solution, inliers] = ransac_point_fit(cor11(:,1:2)', cor22(:,1:2)', 0.01);rmse=0;
% cor1=cor11(inliers,:);
% cor2=cor22(inliers,:);
% fprintf('After RANSAC found %d matches.\n', size(cor1,1));
% button=appendimages(im2,im1,cor2,cor1);

%% FSC
[solution,rmse,cor1,cor2]=FSC(cor11,cor22,'affine',3);
[solution2,rmse2,cor1_2,cor2_2]=FSC(cor11_2,cor22_2,'affine',3);
% 传递坐标数据
aaa = size(cor1,1);
disp(aaa);
bbb = size(cor1_2,1);
disp(bbb);
num2 = aaa+bbb;
disp(num2);
fprintf('After FSC found %d matches.\n', num2);
fprintf('CMR is %f.\n',num2/num1);
fprintf('RMSE is %f.\n', rmse*(aaa/num2)+rmse2*(bbb/num2));
im3=appendimages(im2,im1,cor2,cor1,cor2_2,cor1_2);
cor1=cor1(:,1:6);
cor2=cor2(:,1:6);
cor1_2=cor1_2(:,1:2);
cor2_2=cor2_2(:,1:2);

solution_sum = (solution+solution2)/2;

end