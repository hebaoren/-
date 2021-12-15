function [T] = image_3D_output(image_title, gradient)
%image_3D_output 显示图片的立体信息
% 输入：gradient，任意图像
% 输出：None
    if size(gradient,3)>1
        gradient = rgb2gray(gradient);
    end
    mm = size(gradient,1);
    nn = size(gradient,2);
    max_rc = max(mm,nn);
    if mm>nn
        max_num_cc = blkdiag(gradient,zeros(0,abs(nn-mm)));
    else
        max_num_cc = blkdiag(gradient,zeros(abs(nn-mm),0));
    end
    [ROW,COL] = meshgrid(1:max_rc,1:max_rc);
    figure,meshz(ROW,COL,max_num_cc),title(image_title);
end

