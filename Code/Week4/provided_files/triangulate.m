function Xtrain = triangulate (x1, x2, P1, P2, imsize)
h= [ 2/imsize(1) 0 -1
    0 2/imsize(2) -1
    0 0 1];

x11=euclid(h*homog(x1));
x22=euclid(h*homog(x2));
P11=h*P1;
P22=h*P2;

A = [x11(1)*P11(3,:)-P11(1,:);x11(2)*P11(3,:)-P11(2,:); x22(1)*P22(3,:)-P22(1,:);x22(2)*P22(3,:)-P22(2,:)];

Aprim= [A(:,1) A(:,2) A(:,3)];

[U,D,VT] = svd(Aprim);
UT=transpose(U);
A4prim1=UT*A(:,4);
A4prim=euclid(A4prim1);

for i=1:size(D,2)
    Y(i,1)= (-A4prim(i))/D(i,i);
end

V=transpose(VT);
Xtrain = V*Y;

Xtrain=homog(Xtrain);


end