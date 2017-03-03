function [ outI ] = apply_H( I, H)
% apply_H function that applies transformation H to image I
% returns image outI. 

sizeI = size(I);
sizeH = size(H);
outI = zeros(sizeI(1),sizeI(2),sizeI(3));

% Check dimensions
if ((sizeH(1) ~= 3) || (sizeH(2) ~= 3))
    display('Matrix H has incorrect dimensions')
    return
end

xRange = [0, sizeI(2)-1];
yRange = [0, sizeI(1)-1];

box = [xRange(0); yRange(0); xRange(1); yRange(1)];

% [x,y] = meshgrid(minx:maxx-1,miny:maxy-1);
% 
% pp = p2t(inv(H),[vec(x)';vec(y)']);
% xi=ivec(pp(1,:)',size(x,1));
% yi=ivec(pp(2,:)',size(y,1));
% I2=interp2(0:size(I,2)-1, 0:size(I,1)-1,double(I),xi,yi,meth,NaN);
% 
% if nargout == 3
%     alpha = ~isnan(I2);
% end
% 
% I2 = uint8(I2);




end

