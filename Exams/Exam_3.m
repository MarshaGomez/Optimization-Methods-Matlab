%% 
% Consider the following multi-objective problem:
% 
% $$\left\lbrace \begin{array}{ll}\textrm{minimize} & x_2 \;,-x_1 -2x_2 \\s\ldotp 
% t & x_1 \ge 0\\ & x_2 \ge 0\\ & x_1 +x_2 \le 5\end{array}\right.$$

clear global;
close all;
clear;
clc;

matlab.lang.OnOffSwitchState = 1;
%% 
% *a) Is it a convex problem? Why?*
% 
% We need to analyze the convexity of all the objective function and the feasible 
% region:
%% 
% * Objective Function
% * $f_1 \left(x\right)=x_2$ 
% * $f_2 \left(x\right)={-x}_1 -2x_2$
%% 
% The objective function $f_1 \left(x\right)$ is a linear function, as same 
% as $f_2 \left(x\right)$, which it means both of the objective functions are 
% convex.
%% 
% * Constraints 
% * $g_1 \left(x\right)=x_1 \ge 0$
% * $g_2 \left(x\right)=x_2 \ge 0$
% * $g_3 \left(x\right)=x_1 +x_2 \le 5$
%% 
% The constraints functions $g_1 \left(x\right)\;$and $g_2 \left(x\right)$ are 
% both linear inequalities and the feasible region belong to the non-negative 
% Castersian plane, and it means is convex set.  
% 
% The constraint function $g_3 \left(x\right)$ is also linear inequality constraint 
% and the region is defined by $x_1 +x_2 =5$, $x_1 =0$and $x_2 =0$ makes a triangle 
% set that is also convex.
% 
% Since all the objective functions and the feasible region are convex, we can 
% confirm that this multi-objective optimization problem is convex.

% Definition of the objective functions and the constraints
objective_function_1 = @(x1,x2) x2;
objective_function_2 = @(x1,x2) -x1 - 2.*x2;
constraint_1 = @(x1,x2) -x1;
constraint_2 = @(x1,x2) -x2;
constraint_3 = @(x1,x2) x1 + x2 - 5;

% Creation of the grid with x1 and x2
x1 = linspace(0, 5, 100);
x2 = linspace(0, 5, 100);
[X1,X2] = meshgrid(x1,x2);

% Evaluate the objective functions and the feasible region
F1 = objective_function_1(X1,X2);
F2 = objective_function_2(X1,X2);
C1 = constraint_1(X1,X2);
C2 = constraint_2(X1,X2);
C3 = constraint_3(X1,X2);
%% 
% 

% Create a 3D plot
figure;
surf(X1,X2,F1);
hold on
surf(X1,X2,F2);
surf(X1,X2,C1, 'FaceAlpha', 0.3);
surf(X1,X2,C2, 'FaceAlpha', 0.3);
surf(X1,X2,C3, 'FaceAlpha', 0.3);

xlabel("x1");
ylabel("x2");
zlabel("Objective Function");

legend("F1", "F2", "C1", "C2", "C3");
title("3D Multi-objective Function");

view(3);

hold off
%% 
% 

% Create 2D Plot
figure;
contour(X1,X2,F1, [0 0], 'k', 'LineWidth', 0.5);
hold on

contour(X1,X2,F2, [0 0], 'k', 'LineWidth', 0.5);
contour(X1,X2,C1, [0 0], 'b', 'LineWidth', 0.3);
contour(X1,X2,C2, [0 0], 'g', 'LineWidth', 0.3);
contour(X1,X2,C3, [0 0], 'r', 'LineWidth', 0.3);

legend('F1', 'F2', 'C1', 'C2', 'C3');
xlim([-5 10]);
ylim([-5 10]);
title('2D Multi-objective Optimization Problem');
grid on
hold off
%% 
% *b) Do minima exist? Why?*
% 
% As we see on the previous point, we have a convex optimization problem and 
% the feasible region is in the non-negative quadrant of the cartesian plane that 
% is bounded by a triangular region. The Pareto frontier lies within this feasible 
% region. 
% 
% This given multiobjective optimization problem does have  minima, which correspond 
% to the points on the Pareto frontier within the feasible region defined by the 
% constraints.
% 
% *c) Is the point (0, 2) a weak minimum? Why?*
% 
% First of all, we have
%% 
% * Objective Function
% * $f_1 \left(x\right)=x_2$ 
% * $f_2 \left(x\right)={-x}_1 -2x_2$
%% 
% Let's evaluate both objective functions at point (0,2):

x1 = 0;
x2 = 2;

fprintf("Evaluation at point (%2.2f %2.2f)", x1,x2);
fprintf("F1: %2.2f",objective_function_1(x1,x2));
fprintf("F2: %2.2f",objective_function_2(x1,x2));
%% 
% Now, let's consider neirby point that are within the feasible region and evaluate 
% the objectives again with those points:

epsilon = 0.5;
phi = 0.5;

new_x1 = x1 + epsilon;
new_x2 = x2 - phi;

fprintf("Evaluation at point (%2.2f %2.2f)", new_x1,new_x2);
fprintf("F1: %2.2f",objective_function_1(new_x1,new_x2));
fprintf("F2: %2.2f",objective_function_2(new_x1,new_x2));

new_x1 = x1 - epsilon;
new_x2 = x2 + phi;

fprintf("Evaluation at point (%2.2f %2.2f)", new_x1,new_x2);
fprintf("F1: %2.2f",objective_function_1(new_x1,new_x2));
fprintf("F2: %2.2f",objective_function_2(new_x1,new_x2));
%% 
% It means, exists at least one better objective function values for at least 
% one objective compare to the point (0,2). Therefore, this point is not a weak 
% minima.
% 
% *d) Find all the minima by using the scalarization method.*


%% 
% *e) Find all the weak minima by using the scalarization method.*


%% 
% *f) Find the ideal point.*


%% 
% *g) Apply the goal method with L1*


%% 
% *h) Apply the goal method with L2*


%% 
% *i) Apply the goal method with L-infinity*


%% 
%