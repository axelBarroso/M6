function [ outI ] = apply_H( I, H)
% apply_H function that applies transformation H to image I
% returns image outI. 

sizeI = size(I);
sizeH = size(H);
outI = zeros(sizeI(1),sizeI(2),sizeI(3));
[nrows, ncols, nchan] = size(I);

% Check dimensions
if ((sizeH(1) ~= 3) || (sizeH(2) ~= 3))
    display('Matrix H has incorrect dimensions')
    return
end


xmin = 1;
xmax = ncols;
ymin = 1;
ymax = nrows;

[X,Y] = meshgrid(xmin:xmax, ymin:ymax);
Hncols = xmax - xmin + 1;
Hnrows = ymax - ymin + 1;
Z = ones(Hnrows,Hncols);

XYZs = [X(:) Y(:) Z(:)]';

Hi = inv(H); 
HiXYZs = Hi * XYZs;
HX = reshape(HiXYZs(1,:), Hnrows, Hncols);
HY = reshape(HiXYZs(2,:), Hnrows, Hncols);
HZ = reshape(HiXYZs(3,:), Hnrows, Hncols);
HX = HX ./ HZ;
HY = HY ./ HZ;

nchan= size(I,3);

for c=1:nchan,
    outI(:,:,c) = interp2(double(I(:,:,c)), HX, HY, 'linear', 0);
end

figure; imshow(uint8(I)); title('original image');
figure; imshow(uint8(outI)); title('transformed image');

end