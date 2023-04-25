%% 
% *Multiobjective optimization. Goal Method.*
% 
% $$\left\lbrace \begin{array}{ll}\textrm{minimize} & \left(x_1 +2x_2 -3x_3 
% \right),\left(-x_1 -x_2 -x_3 \right),\left(-4x_1 -2x_2 +x_3 \right)\\\textrm{subject}\;\textrm{to} 
% & x_1 +x_2 +x_3 \le 10\\\; & x_3 \le 5\\\; & x_1 \ge 0\\\; & x_2 \ge 0\\\; & 
% x_3 \ge 0\end{array}\right.$$
% 
% a) Find the ideal point. 

close all;
clear;
clc;
matlab.lang.OnOffSwitchState = 1;


C = [ 1  2 -3
     -1 -1 -1
     -4 -2  1];

A = [ 1  1  1
      0  0  1
     -1  0  0
      0 -1  0
      0  0 -1];

b = [10
      5
      0
      0
      0];

n = size(C,2); 
p = size(C,1);
m = size(A,1);

% Find the ideal point

options = optimset('Display', 'off');
z = zeros(p,1);

for i = 1:p
    [~,z(i)] = linprog(C(i,:), A, b, [],[],[],[],options);
end

%% 
% b) Apply the goal method with s = 1. 

c = [zeros(n,1)
     ones(p,1)];

P = [C -eye(p)
    -C -eye(p)
     A zeros(m,p)];

q = [z
    -z
     b];

solution_1 = linprog(c, P,q, [],[],[],[],options);
y_1 = solution_1(1:n);
%% 
% c) Apply the goal method with s = 2. 

H = C'*C;
f = C'*z;

solution_2 = quadprog(H,f,A,b, [], [], [],[],[],options);
y_2 = solution_2(1:n);
%% 
% d) Apply the goal method with s = +∞. Is the found point a minimum?

c = [zeros(n,1)
     1];

P = [C -ones(p,1)
    -C -ones(p,1)
     A zeros(m,1)];

q = [z
    -z
     b];

solution_inf = linprog(c,P,q,[],[],[],[],options);
y_inf = solution_inf(1:n);
%%
c = [zeros(n,1) ; -ones(p,1)] ;
P = [C eye(p); 
     A zeros(m,p) ; 
     zeros(n,n) -eye(p)] ;
q = [C*y_inf ; b ; zeros(p,1)] ;
[~,v_minimum] = linprog(c,P,q,[],[],[],[],[],options) 
%% 
% 
% 
% 
% 
%