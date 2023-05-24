%% 
% Consider the following optimization problem:
% 
% $$\left\lbrace \begin{array}{ll}\mathrm{minimize} & 2x_1^2 +x_1 x_3 +x_2^2 
% +2x_2 x_3 +3x_3^2 +x_3 x_4 +2x_4^2 -5x_1 -4x_3 +3x_4 \\x\;\epsilon \;\Re^4  
% & \;\end{array}\right.$$

close all;
clear;
clc;

matlab.lang.OnOffSwitchState = 1;
%% 
% *a) Do global optimal solutions exist? Why?*


%% 
% *b) Is it a convex problem? Why?*
% 
% We need to analyze the covexity of each term and the objective function as 
% a whole. A function is convex if the Hessian Matrix is positive semi-definite. 
% 
% In this case the objective function is:
% 
% $$f\left(x\right)\;=\;2x_1^2 +x_1 x_3 +x_2^2 +2x_2 x_3 +3x_3^2 +x_3 x_4 +2x_4^2 
% -5x_1 -4x_3 +3x_4$$
% 
% To create the Hessian Matrix in a quick form we have:
% 
% General Form = $a_1 x_1^2 +a_2 x_2^2 +a_3 x_3^2 \;\;{+a_4 x}_4^2 +a_5 x_1 
% x_2 \;+a_6 x_1 x_3 \;+a_7 x_1 x_4 +a_8 x_2 x_3 \;+a_9 x_2 x_4 +a_{10} x_3 x_4$
% 
% Hessian = $\left\lbrack \begin{array}{cccc}2a_1  & a_5  & a_6  & a_7 \\a_5  
% & 2a_2  & a_8  & a_9 \\a_6  & a_8  & 2a_3  & a_{10} \\a_7  & a_9  & a_{10}  
% & 2a_4 \end{array}\right\rbrack$
% 
% To check if the Hessian Matrix is semi-definite, we need to compute the eigenvalues 
% of the matrix.

H = [4 0 1 0
     0 2 2 0
     1 2 6 1
     0 0 1 4];

fprintf('Eigen Values Objective Function:');
display(eig(H))
%% 
% Since all the eigenvalues are positive and not zero, we can conclude that 
% Hessian Matrix is positive definite, therefore is a strongly convex optimization 
% problem.
% 
% *c) Find the global minimum by using the gradient method with exact line search.*
% 
% Starting from the point $\left(0,0,0,0\right)\;\;\mathrm{with}\;\left\|\nabla 
% f\left(x\right)\right\|<{10}^{-6}$ as stopping criterion. How many iterations 
% are needed?


%% 
% *d) Find the global minimum by using the conjugate gradient method* 
% 
% Starting from the point (0, 0, 0, 0). How many iterations are needed? Write 
% the point found by the method at any iteration.


%% 
% *e) Find the global minimum by using the Newton method* 
% 
% Starting from the point (0, 0, 0, 0). How many iterations are needed?