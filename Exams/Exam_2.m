%% 
% Consider the following optimization problem:
% 
% $$\left\lbrace \begin{array}{ll}\textrm{minimize} & 2x_1^2 +x_1 x_3 +x_2^2 
% +2x_2 x_3 +3x_3^2 +x_3 x_4 +2x_4^2 -5x_1 -4x_3 +3x_4 \\x\;\epsilon \;\Re^4  
% & \;\end{array}\right.$$

close all;
clear;
clc;

matlab.lang.OnOffSwitchState = 1;
%% 
% *a) Do global optimal solutions exist? Why?*
% 
% To determine convexity, we need to examine the Hessian Matrix, the second 
% partial derivates of the objective function with respect to the variables. On 
% the next 1.2 we can check that our Hessian Matrix is strongly convex and we 
% know that in any convex optimization problem the feasible region is a convex 
% set. 
% 
% Using the Theorem weierstrass that says if the objective funciton is continuous 
% and the feasible region is closed and bounded, then at least a global optimum 
% exists. 
% 
% On this problem using Collory 3 were the objective function is continuous 
% and coercive and the feasible region is closed, then at least a global minimum 
% exists. Coercivity is$\lim_{\left\|x\right\|\to \infty } f\left(x\right)=+\infty$

% In this case we use 1 variable, because is an approximation calculation
objective_function = @(x) 2*x^2 + x*x + x^2 + 2*x*x + 3*x^2 + x*x + 2*x^2 - 5*x - 4*x + 3*x;
syms x

limit = limit(objective_function(x), x, Inf)
%% 
% Since the objective function is coercive, there exists a global optimum.
% 
% *b) Is it a convex problem? Why?*
% 
% We need to analyze the convexity of each term and the objective function as 
% a whole. A function is convex if the Hessian Matrix is positive semi-definite. 
% 
% In this case the objective function is:
% 
% $$f\left(x\right)\;=\;2x_1^2 +x_1 x_3 +x_2^2 +2x_2 x_3 +3x_3^2 +x_3 x_4 +2x_4^2 
% -5x_1 -4x_3 +3x_4$$
% 
% To create the Hessian Matrix in a quick form we have:
% 
% $$\textrm{General}\;\textrm{Form}=a_1 x_1^2 +a_2 x_2^2 +a_3 x_3^2 \;\;+a_4 
% x_4^2 \;\;+a_5 x_1 x_2 \;+a_6 x_1 x_3 \;+a_7 x_1 x_4 +a_8 x_2 x_3 \;+a_9 x_2 
% x_4 +a_{10} x_3 x_4$$
% 
% $$\textrm{Hessian}=\left\lbrack \begin{array}{cccc}2a_1  & a_5  & a_6  & a_7 
% \\a_5  & 2a_2  & a_8  & a_9 \\a_6  & a_8  & 2a_3  & a_{10} \\a_7  & a_9  & a_{10}  
% & 2a_4 \end{array}\right\rbrack$$
% 
% To check if the Hessian Matrix is semi-definite, we need to compute the eigenvalues 
% of the matrix.

Q = [4 0 1 0
     0 2 2 0
     1 2 6 1
     0 0 1 4];

fprintf('Eigen Values Objective Function:');
display(eig(Q))
fprintf('Determinant: %i', det(Q));
%% 
% Since all the eigenvalues are positive and not zero, also the Determinant 
% of the Hessian matrix positive definite, we can conclude that Hessian Matrix 
% is positive definite, therefore is a strongly convex optimization problem.
% 
% *c) Find the global minimum by using the gradient method with exact line search.*
% 
% Starting from the point $\left(0,0,0,0\right)\;\;\textrm{with}\;\left\|\nabla 
% f\left(x\right)\right\|<{10}^{-6}$ as stopping criterion. How many iterations 
% are needed?
% 
% The optimization problem is unconstraint and the objective function is a quadratic, 
% where Q is positive definite, as follows:
% 
% $$\begin{array}{l}f\left(x\right)=\frac{1}{2}x^T Q\;x+q^T x\\g\left(x\right)=Q\;x+q;\end{array}$$
% 
% The steepest descent direction:
%%
% 
%   % Gradient method
%   Choose x0 tolerance=ie-6 i=0
%   while (norm(g) > tolerance)
%      [step direction] d = -g;
%      [step size] alpha = -(g'*d)/d'Qd;
%      x = x + alpha*d;
%      i = i + 1;
%   end
%

q = [-5
      0
     -4
      3];

x0 = [0
      0
      0
      0];

tolerance = 1e-6;

% Starting Point
x = x0;
i = 0;

g = Q*x + q;

fprintf('Gradient Method with exact line search');
fprintf('i \t f(x) \t \t \t ||grad f(x)|| \n');

while (norm(g)>tolerance)
    v = (1/2)*x'*Q*x + q'*x;
    g = Q*x + q;
   
    fprintf('%i \t %2.10f \t \t %2.10f \n',i,v,norm(g));
    
    % Direction
    d = -g;
    % Step size
    alpha = -(g'*d)/(d'*Q*d);
    % New point
    x = x + alpha*d;
    % Next iteration
    i = i+1;
end
%% 
% With Gradient Method using exact line search we need 48 iterations to arrive 
% to the best solution.
% 
% *d) Find the global minimum by using the conjugate gradient method* 
% 
% Starting from the point (0, 0, 0, 0). How many iterations are needed? Write 
% the point found by the method at any iteration. With this method, the search 
% direction involves the gradient computed at the current iteration and the direction 
% computed at the precious iteration. 
% 
% $$d^k =\left\lbrace \begin{array}{ll}-g^0  & \textrm{if}\;k=0\\-g^k +\beta_k 
% d^{k-1}  & \textrm{if}\;k\ge 1\end{array}\right.$$
% 
% The optimization problem is unconstraint and the objective function is a quadratic, 
% where Q is positive definite, as follows:
% 
% $$\begin{array}{l}f\left(x\right)=\frac{1}{2}x^T Q\;x+q^T x\\g\left(x\right)=Q\;x+q;\end{array}$$
% 
% Conjugate gradient method for quadratic functions with linear search:
%%
% 
%   % Conjugate method
%   Choose x0 i=0
%   while g > 0
%      if i = 0
%          d = -g;
%      else 
%          beta = norm(g)^2 / norm(g_prev)^2;
%          d = -g + beta*d_prev;
%      end
%      alpha = norm(g)^2 / (d'Q*d);
%      x = x + alpha*d;
%      g = Qx + q;
%      i = i + 1;
%   end
%

clear x v g

x0 = [0
      0
      0
      0];

x = x0;
i = 0;

v = (1/2)*x'*Q*x + q'*x;
g = Q*x + q;

fprintf("\n \n Conjugate Gradient with exact line search \n");
fprintf("i \t f(x) \t \t \t ||grad f(x)|| ");

while (norm(g)>tolerance)
    fprintf('%i \t %2.10f \t \t %2.10f \n', i, v, norm(g));
    if i == 0
        % Direction on the first iteration
        d = -g;
    else 
        beta = (norm(g)^2) / (norm(g_prev)^2);
        % Directio
        d = -g + beta*d_prev;
    end

    % Step Size
    alpha = (norm(g)^2) / (d'*Q*d);
    x = x + alpha*d;
   
    d_prev = d;
    g_prev = g;

    v = (1/2)*x'*Q*x + q'*x;
    g = Q*x + q;

    i = i + 1;
end
%% 
% On this case, we can easily see that the convergence rate is quicklier than 
% the previous method, with Conjugate Method using exact line search we use just 
% 3 iterations to found the minimum.
% 
% *e) Find the global minimum by using the Newton method* 
% 
% Starting from the point (0, 0, 0, 0). How many iterations are needed?
% 
% We can conclude an objective function is strongly convex because the Hessian 
% of the function is positive definite, and we get the global convergence because 
% the direction is a descent direction:
% 
% $$\nabla f{\left(x^k \right)}^T d^k =-\nabla f{\left(x^k \right)}^T {\left\lbrack 
% \nabla^2 f\left(x^k \right)\right\rbrack }^{-1} \nabla \left(f\left(x^k \right)\right)<0$$
% 
% Newton method basic version:
%%
% 
%   % Newton method
%   Choose x0 i=0
%   while (g ~= 0)
%       [search direction] d = -H\g;
%       x = x + d 
%       i = i + 1;
%   end
%
%% 
% 

x0 = [0
      0
      0
      0];

i = 0;
x = x0;

v = (1/2)*x'*Q*x + q'*x;
g = Q*x + q;

fprintf("\n \n Newton method basic version \n");
fprintf("i \t f(x) \t \t \t ||grad f(x)|| \n");
while (norm(g) > tolerance)
    fprintf("%i \t %2.10f \t \t %2.10f \n", i, v, norm(g));
    d = -Q\g;
    x = x + d;
    i = i + 1;

    v = (1/2)*x'*Q*x + q'*x;
    g = Q*x + q;
end