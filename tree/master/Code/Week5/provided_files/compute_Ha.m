function [ Ha ] = compute_Ha(Pproj, Hp, v1, v2, v3 )
% ToDo: compute the matrix Ha that 
%       upgrades the projective reconstruction to an affine reconstruction
% Use the following vanishing points given by three pair of orthogonal lines
% and assume that the skew factor is zero and that pixels are square

    % From lecture 9, page 35/40
    A = [v1(1)*v2(1) v1(1)*v2(2)+v1(2)*v2(1) v1(1)*v2(3)+v1(3)*v2(1) v1(2)*v2(2) v1(2)*v2(3)+v1(3)*v2(2) v1(3)*v2(3); ...
        v1(1)*v3(1) v1(1)*v3(2)+v1(2)*v3(1) v1(1)*v3(3)+v1(3)*v3(1) v1(2)*v3(2) v1(2)*v3(3)+v1(3)*v3(2) v1(3)*v3(3); ...
        v2(1)*v3(1) v2(1)*v3(2)+v2(2)*v3(1) v2(1)*v3(3)+v2(3)*v3(1) v2(2)*v3(2) v2(2)*v3(3)+v2(3)*v3(2) v2(3)*v3(3); ...
        0 1 0 0 0 0; ...
        1 0 0 -1 0 0];

    s = null(A);

    s = s(:,3);
    wmetric = [s(1) s(2) s(3); ...
        s(2) s(4) s(5); ...
        s(3) s(5) s(6)];

    p1 = Pproj(1:3, :);
    p2 = Pproj(4:end, :);

    pm = p2 * inv(Hp);
    pme = pm(1:3, 1:3);

    % use the Cholesky decomposition of S to compute an upper
    % triangular matrix K such that S = KK^T :
    K  = chol(inv(pme' * wmetric * pme));

    % the matrix K is a possible matrix A that can be used to metrically
    % rectify the image.
    Ha = [inv(K), [0, 0, 0]'; 0, 0, 0, 1];

end

