%% 
% *Convexity* 
% 
% Prove that if C is convex, then for any $x^1 ,\ldotp \ldotp \ldotp ,x^k \in 
% C$and $\alpha_1 ,\ldotp \ldotp \ldotp ,\alpha_k \in \left(0,1\right)\;$subject 
% to$\sum_{i=1}^k \alpha_i =1\;$one has $\sum_{i=1}^k \alpha_i x^i \in C$

close all;
clear;
clc;

% Define the problem statement
% We want to prove that if C is convex, then for any x_1, ..., x_k ∈ C and α_1, ..., α_k ∈ (0, 1),
% the point y = α_1*x_1 + ... + α_k*x_k also belongs to C

% Generate some random points for C
C = rand(10,2);

% Generate some random weights for α_1, ..., α_k
alpha = rand(1,10);

% Calculate the point y = α_1*x_1 + ... + α_k*x_k
y = alpha * C;

% Check if y belongs to C by checking if all the line segments between x_i and y lie within C
for i = 1:size(C,1)
    % Calculate the line segment between x_i and y
    segment = [C(i,:); y];
    
    % Check if the line segment lies within C using the "inpolygon" function
    if ~inpolygon(segment(:,1), segment(:,2), C(:,1), C(:,2))
        % If the line segment does not lie within C, print an error message and break out of the loop
        fprintf('Error: y does not belong to C\n');
        return
    end
end

% If we get to this point, y belongs to C for the given x_i and α_i
fprintf('y belongs to C');