%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 1: Image rectification


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Applying image transformations

% ToDo: create the function  "apply_H" that gets as input a homography and
% an image and returns the image transformed by the homography.
% At some point you will need to interpolate the image values at some points,
% you may use the Matlab function "interp2" for that.


%% 1.1. Similarities
I=imread('Data/0000_s.png'); % we have to be in the proper folder
figure(1); imshow(I); title('Original Image');

% ToDo: generate a matrix H which produces a similarity transformation
% Matrix H with translation
H = [1 0 100; 
     0 1 150;
     0 0  1 ];

I2 = apply_H(I, H, 'keepOriginalPositions');
figure(2); imshow(uint8(I2)); title('Planar transformations: Translation');

% Matrix H with rotation of 20º
H = [cosd(20) -sind(20) 0; 
     sind(20) cosd(20) 0;
     0 0 1];

I2 = apply_H(I, H);
figure(3); imshow(uint8(I2)); title('Planar transformations: Rotation');

% Matrix H with translation & rotation
H = [cosd(20) -sind(20) 100; 
     sind(20) cosd(20) 150;
     0 0 1];

I2 = apply_H(I, H, 'keepOriginalPositions');
figure(4); imshow(uint8(I2)); title('Planar transformations: Translation & Rotation');

% Matrix H with scale
H = [2 0 0; 
     0 2 0;
     0 0 1];

I2 = apply_H(I, H);
figure(5); imshow(uint8(I2)); title('Planar transformations: Scale');


%% 1.2. Affinities

% ToDo: generate a matrix H which produces an affine transformation
H = [1.15  -0.3  160; 
     0.3   0.65  -65;
     0      0     1];

I2 = apply_H(I, H);
figure(6); imshow(uint8(I2)); title('Planar transformations: Affine')

% ToDo: decompose the affinity in four transformations: two
% rotations, a scale, and a translation
A = H(1:size(H,1)-1,1:size(H,2)-1);
[U,S,V] = svd(A); %performs a singular value decomposition of matrix A, such that A = U*S*V'.
R_theta = U*V';
R_phi = V';
D = S;
t = [H(1,size(H,1)); H(2,size(H,1))];

% ToDo: verify that the product of the four previous transformations
% produces the same matrix H as above

% Compute A, following Lecture 1, Page 43 of 48
A = R_theta*R_phi'*D*R_phi;
% Add the translation terms to H
H_new = [A t; 0 0 1];
I2 = apply_H(I, H_new);
figure(7); imshow(uint8(I2)); title('Planar transformations: Affine (From SVD decomposition)')

% ToDo: verify that the proper sequence of the four previous
% transformations over the image I produces the same image I2 as before
Identity = [1 0; 0 1];
H_R_phi = [R_phi [0; 0]; 0 0 1];
H_D = [D [0; 0]; 0 0 1];
H_R_theta = [R_theta [0; 0]; 0 0 1];
H_t = [Identity t; 0 0 1];

I_aux = apply_H(I, H_R_theta);
I_aux = apply_H(I_aux, H_R_phi');
I_aux = apply_H(I_aux, H_D);
I_aux = apply_H(I_aux, H_R_phi);
I_aux = apply_H(I_aux, H_t);
figure(8); imshow(uint8(I_aux)); title('Planar transformations: Affine (Apply matrices independently)')

%% 1.3 Projective transformations (homographies)
% ToDo: generate a matrix H which produces a projective transformation
I2 = apply_H(I, H);
figure(9); imshow(uint8(I2));title('Projective transformations: Not DONE !! ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Affine Rectification

% choose the image points
I = imread('Data/0000_s.png');
A = load('Data/0000_s_info_lines.txt');

% indices of lines
i = 424;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 240;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 712;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 565;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';
i = 227;
p9 = [A(i,1) A(i,2) 1]';
p10 = [A(i,3) A(i,4) 1]';
i = 576;
p11 = [A(i,1) A(i,2) 1]';
p12 = [A(i,3) A(i,4) 1]';
i = 534;
p13 = [A(i,1) A(i,2) 1]';
p14 = [A(i,3) A(i,4) 1]';

% ToDo: compute the lines l1, l2, l3, l4, that pass through the different pairs of points
l1 = computeLine( p1, p2);
l2 = computeLine( p3, p4);
l3 = computeLine( p5, p6);
l4 = computeLine( p7, p8);

% show the chosen lines in the image
figure(10);imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');
% plot(t, -(l5(1)*t + l5(3)) / l5(2), 'y');
% plot(t, -(l6(1)*t + l6(3)) / l6(2), 'y');

% ToDo: compute the homography that affinely rectifies the image
% Compute vanishing points
vp_12 = cross(l1, l2);
vp_12 = vp_12 / vp_12(3);
vp_34 = cross(l3, l4);
vp_34 = vp_34 / vp_34(3);
l_inf = cross(vp_12,vp_34);
l_inf = l_inf / l_inf(3);

H_ap = [1 0 0; 0 1 0; l_inf];
I_ap = apply_H(I, H_ap);
figure(11); imshow(uint8(I_ap)); title('Affine rectification via the vanishing line')

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4
new_p1  = computeTransformedPoints(H_ap, p1);
new_p2  = computeTransformedPoints(H_ap, p2);
new_p3  = computeTransformedPoints(H_ap, p3);
new_p4  = computeTransformedPoints(H_ap, p4);
new_p5  = computeTransformedPoints(H_ap, p5);
new_p6  = computeTransformedPoints(H_ap, p6);
new_p7  = computeTransformedPoints(H_ap, p7);
new_p8  = computeTransformedPoints(H_ap, p8);
new_p9  = computeTransformedPoints(H_ap, p9);
new_p10 = computeTransformedPoints(H_ap, p10);
new_p11 = computeTransformedPoints(H_ap, p11);
new_p12 = computeTransformedPoints(H_ap, p12);
new_p13 = computeTransformedPoints(H_ap, p13);
new_p14 = computeTransformedPoints(H_ap, p14);

lr1 = computeLine( new_p1, new_p2);
lr2 = computeLine( new_p3, new_p4);
lr3 = computeLine( new_p5, new_p6);
lr4 = computeLine( new_p7, new_p8);
lr5 = computeLine( new_p9, new_p12);
lr6 = computeLine( new_p10, new_p13);

% show the transformed lines in the transformed image
figure(12);imshow(uint8(I_ap)); title('Affine rectification. Rectificated lines')
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'y');

% ToDo: to evaluate the results, compute the angle between the different pair 
% of lines before and after the image transformation

% get normalized version of lines (two coordinates)
norml1 = [l1(1)/l1(3), l1(2)/l1(3)];
norml2 = [l2(1)/l2(3), l2(2)/l2(3)];
norml3 = [l3(1)/l1(3), l3(2)/l3(3)];
norml4 = [l4(1)/l4(3), l4(2)/l4(3)];

normlr1 = [lr1(1)/lr1(3), lr1(2)/lr1(3)];
normlr2 = [lr2(1)/lr2(3), lr2(2)/lr2(3)];
normlr3 = [lr3(1)/lr1(3), lr3(2)/lr3(3)];
normlr4 = [lr4(1)/lr4(3), lr4(2)/lr4(3)];

angle13 = acosd(dot(norml1,norml3)/(norm(norml1)*norm(norml3)));
angler13 = acosd(dot(normlr1,normlr3)/(norm(normlr1)*norm(normlr3)));

angle14 = acosd(dot(norml1,norml4)/(norm(norml1)*norm(norml4)));
angler14 = acosd(dot(normlr1,normlr4)/(norm(normlr1)*norm(normlr4)));

angle23= acosd(dot(norml2,norml3)/(norm(norml2)*norm(norml3)));
angler23 = acosd(dot(normlr2,normlr3)/(norm(normlr2)*norm(normlr3)));

angle24= acosd(dot(norml2,norml4)/(norm(norml2)*norm(norml4)));
angler24 = acosd(dot(normlr2,normlr4)/(norm(normlr2)*norm(normlr4)));

disp(['Upper left corner before transformation: ' , num2str(angle13)]);
disp(['Upper left corner after transformation: ' , num2str(angler13)]);
disp(' ');

disp(['Upper right corner before transformation: ' , num2str(angle14)]);
disp(['Upper right corner after transformation: ' , num2str(angler14)]);
disp(' ');

disp(['Lower left corner before transformation: ' , num2str(angle23)]);
disp(['Lower left corner after transformation: ' , num2str(angler23)]);
disp(' ');

disp(['Lower right corner before transformation: ' , num2str(angle24)]);
disp(['Lower right corner after transformation: ' , num2str(angler24)]);
disp(' ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Metric Rectification

%% 3.1 Metric rectification after the affine rectification (stratified solution)

% ToDo: Metric rectification (after the affine rectification) using two non-parallel orthogonal line pairs
%       As evaluation method you can display the images (before and after
%       the metric rectification) with the chosen lines printed on it.
%       Compute also the angles between the pair of lines before and after
%       rectification.
lr5 = computeLine( new_p9, new_p12);
lr6 = computeLine( new_p10, new_p13);
lr7 = computeLine( new_p2, new_p5);
lr8 = computeLine( new_p1, new_p8);

% normalize orthogonal lines so that they have 2 elements
lr5 = lr5 / lr5(3);
lr6 = lr6 / lr6(3);
lr7 = lr7 / lr7(3);
lr8 = lr8 / lr8(3);

figure(10);imshow(uint8(I_ap));
hold on;
t=1:0.1:1000;
plot(t, -(lr5(1)*t + lr5(3)) / lr5(2), 'y');
plot(t, -(lr6(1)*t + lr6(3)) / lr6(2), 'y');
plot(t, -(lr7(1)*t + lr7(3)) / lr7(2), 'y');
plot(t, -(lr8(1)*t + lr8(3)) / lr8(2), 'y');

n = [lr5(1)*lr6(1), lr5(1)*lr6(2) + lr5(2)*lr6(1), lr5(2)*lr6(2) ; ...
     lr7(1)*lr8(1), lr7(1)*lr8(2) + lr7(2)*lr8(1), lr7(2)*lr8(2)];
s = null(n);
S = [s(1) s(2);s(2) s(3)];
% Cholesky decomposition
K = chol(S,'upper');

H_sa = [K [0; 0]; 0 0 1];
H_sa = inv(H_sa);
I2 = apply_H(I_ap, H_sa);
figure(13); imshow(uint8(I2)); title('Metric rectification via orthogonal lines.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. OPTIONAL: Metric Rectification in a single step
% Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)

%% 5. OPTIONAL: Affine Rectification of the left facade of image 0000

%% 6. OPTIONAL: Metric Rectification of the left facade of image 0000

%% 7. OPTIONAL: Affine Rectification of the left facade of image 0001

%% 8. OPTIONAL: Metric Rectification of the left facade of image 0001