%% 
% *Gradient Method. Exact Line Search*
% 
% Implement in MATLAB the gradient method for solving the problem
% 
% $$\left\lbrace \begin{array}{ll}\textrm{minimize} & \frac{1}{2}x^T \textrm{Qx}+c^T 
% x\\x\;\epsilon \;R^n  & \;\end{array}\right.$$
% 
% Solve the Problem
% 
% $$\left\lbrace \begin{array}{ll}\textrm{minimize} & 3x_1^2 +3x_2^2 +3x_3^2 
% +3x_4^2 -4x_1 x_3 -4x_2 x_4 +x_1 -x_2 +2x_3 -3x_4 \\x\;\epsilon \;R^4  & \;\end{array}\right.$$

% Define the objective to be minimize 
function f = objective_function(x):
    % Quadratic function
    f = (1/2).*x'.*Q.*x + c'.*x;
end
%%
% Define the gradient function of the objective function
function g = gradient_function(x):
    
end
%%
% Define the line search function
%%
% Define the gradient method with exact line search