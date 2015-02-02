% Checking possibilities of taking good set of inliers in VP classification
% Number of lines for orthogonal VPs: 1-x, 2-y, 3-z, 4-outliers
n = [10 10 10 10];

% Number of possible combinations taking 2+1 set
combs_2_1 = nchoosek(n(1),2) * (n(2)+n(3)) + ...
            nchoosek(n(2),2) * (n(3)+n(1)) + ...
            nchoosek(n(3),2) * (n(1)+n(2));
        
combs_1_1_1 = prod(n(1:3));

% All possible combinations
all_combs = nchoosek(sum(n),3);

% Probability of taking a set of inliers only
wn = (combs_2_1+combs_1_1_1) / all_combs; % With 1-1-1 combinations
wn_= combs_2_1 / all_combs; % Without 1-1-1 combinations

% Number of necessary steps to assure p success probability
p = 0.9;
k = log(1-p)/log(1-wn);
k_= log(1-p)/log(1-wn_);

fprintf('Number of 2-1 combinations:\t%d\n',combs_2_1);
fprintf('Number of 1-1-1 combinations:\t%d\n',combs_1_1_1);
fprintf('All possible combinations:\t%d\n',all_combs);
fprintf('Probability of taking set of inliers: %.2f\n', wn);
fprintf('Number of steps for %.2f success probability (with 1-1-1): %d\n', p, ceil(k));
fprintf('Number of steps for %.2f success probability (no   1-1-1): %d\n', p, ceil(k_));

% Computational cost
fprintf('#Cost (with 1-1-1): %d\n', ceil(k)*4);
fprintf('#Cost (no   1-1-1): %d\n', ceil(k_)*3);