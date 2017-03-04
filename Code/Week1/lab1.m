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

% ToDo: generate a matrix H which produces a similarity transformation
% Matrix H with translation
H = [1 0 100; 
     0 1 150;
     0 0 1];

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));

% Matrix H with rotation of 20º
H = [cosd(20) -sind(20) 0; 
     sind(20) cosd(20) 0;
     0 0 1];

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));

% Matrix H with translation & rotation
H = [cosd(20) -sind(20) 100; 
     sind(20) cosd(20) 150;
     0 0 1];

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));

% Matrix H with scale
H = [2 0 0; 
     0 2 0;
     0 0 1];

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));


%% 1.2. Affinities

% ToDo: generate a matrix H which produces an affine transformation
H = [1.15 -0.3 160; 
     0.3 0.65 -65;
     0 0 1];

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));

% ToDo: decompose the affinity in four transformations: two
% rotations, a scale, and a translation
[U,S,V] = svd(H(1:size(H,1)-1,1:size(H,2)-1)); %performs a singular value decomposition of matrix A, such that A = U*S*V'.
l = [0,0,1];
padding = [0,0]';
rotation_1 = [U*V' padding; l];
rotation_2 = [V' padding; l];
scaling = [S padding; l];
translation = [1 0 H(1,size(H,1)); 0 1 H(2,size(H,1)); l];


% ToDo: verify that the product of the four previous transformations
% produces the same matrix H as above

% Compute A, following Lecture 1, Page 43 of 48
A = rotation_1*rotation_2'*scaling*rotation_2;
% Add the translation terms to H
H_new = translation*A;
I2 = apply_H(I, H_new);
figure; imshow(uint8(I2));


% ToDo: verify that the proper sequence of the four previous
% transformations over the image I produces the same image I2 as before
I_aux = apply_H(I, rotation_2);
I_aux = apply_H(I_aux, scaling);
I_aux = apply_H(I_aux, rotation_2');
I_aux = apply_H(I_aux, rotation_1);
I_aux = apply_H(I_aux, translation);
figure; imshow(uint8(I_aux));


%% 1.3 Projective transformations (homographies)

% ToDo: generate a matrix H which produces a projective transformation

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Affine Rectification


% choose the image points
I=imread('Data/0000_s.png');
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

% ToDo: compute the lines l1, l2, l3, l4, that pass through the different pairs of points


% show the chosen lines in the image
figure;imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');

% ToDo: compute the homography that affinely rectifies the image

I2 = apply_H(I, H);
figure; imshow(uint8(I2));

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4

% show the transformed lines in the transformed image
figure;imshow(uint8(I2));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'y');

% ToDo: to evaluate the results, compute the angle between the different pair 
% of lines before and after the image transformation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Metric Rectification

%% 3.1 Metric rectification after the affine rectification (stratified solution)

% ToDo: Metric rectification (after the affine rectification) using two non-parallel orthogonal line pairs
%       As evaluation method you can display the images (before and after
%       the metric rectification) with the chosen lines printed on it.
%       Compute also the angles between the pair of lines before and after
%       rectification.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. OPTIONAL: Metric Rectification in a single step
% Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)

%% 5. OPTIONAL: Affine Rectification of the left facade of image 0000

%% 6. OPTIONAL: Metric Rectification of the left facade of image 0000

%% 7. OPTIONAL: Affine Rectification of the left facade of image 0001

%% 8. OPTIONAL: Metric Rectification of the left facade of image 0001