function [F, inliers] = ransac_fundamental_matrix( p1, p2, th )
%RANSAC_FUNDAMENTAL_MATRIX 


[Ncoords, Npoints] = size(p1);

% ransac
it = 0;
best_inliers = [];
% probability that at least one random sample set is free of outliers
p = 0.999; 
max_it = 1000;
while it < max_it
    
    points = randomsample(Npoints, 4);
    F_estimation = fundamental_matrix(p1(:,points), p2(:,points));    
    inliers = compute_inliers(F_estimation, p1, p2, th);
    
    % test if it is the best model so far
    if length(inliers) > length(best_inliers)
        best_inliers = inliers;
    end    
    
    % update estimate of max_it (the number of trials) to ensure we pick, 
    % with probability p, an initial data set with no outliers
    fracinliers =  length(inliers)/Npoints;
    pNoOutliers = 1 -  fracinliers^4;
    pNoOutliers = max(eps, pNoOutliers);  % avoid division by -Inf
    pNoOutliers = min(1-eps, pNoOutliers);% avoid division by 0
    max_it = log(1-p)/log(pNoOutliers);
    
    it = it + 1;
end

% compute F from all the inliers
F = fundamental_matrix(p1(:,best_inliers), p2(:,best_inliers));  
inliers = best_inliers;
end


function idx_inliers = compute_inliers(F, p1, p2, th)
    
    error = [];
    for i = 1 : size(p1, 2)
        % Compute Sampson error 
        error = [error compute_Sampson_error(F, p1, p2)];
    end
    
    idx_inliers = find(error < th.^2);
end

function error = compute_Sampson_error(F, p1, p2)
% Sampson error (1st order approx. of the geometric distance)




end

