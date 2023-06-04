%% 
% 
% 
% 

Q = [4 0 0
     0 4 0
     0 0 4];

eig(Q)
%%

syms x1 x2 x3
objective_function = (x1 - 7)^2 + (x2 - 6)^2 + (x3 - 3)^2;

d_x1 = diff(objective_function, x1)
d_x2 = diff(objective_function, x2)
d_x3 = diff(objective_function, x3)

dd_x1 = diff(d_x1, x1)
dd_x2 = diff(d_x2, x2)
dd_x3 = diff(d_x3, x3)
dd_x1x2 = diff(dd_x1, x2)
dd_x1x3 = diff(dd_x1, x3)
dd_x2x3 = diff(dd_x2, x3)


H = [dd_x1   dd_x1x2 dd_x1x3
     dd_x1x2 dd_x2   dd_x2x3
     dd_x1x3 dd_x2x3 dd_x2]

constraint_active= x1+x2+x3-10;
g_x1 = diff(constraint_active, x1)
g_x2 = diff(constraint_active, x2)
g_x3 = diff(constraint_active, x3)

%% 
% $$\begin{array}{l}0\le x_1 \le 10\\0\le x_2 \le 10\\0\le x_3 \le 10\end{array}$$
% 
% Closed and bounded. Objective function continuous. Strongly convex



f1 = @(x1,x2,x3) (x1 - 7)^2 + (x2 - 6)^2 + (x3 - 3)^2;
c1=@(x1,x2,x3) x1+x2+x3-10; 
c2=@(x1,x2,x3) -x1;
c3=@(x1,x2,x3) -x2;
c4=@(x1,x2,x3) -x3;

x1 = linspace(0,10, 100);
x2 = linspace(0,10, 100);
x3 = linspace(0,10, 100);

[X1,X2,X3] = meshgrid(x1,x2,x3);

C1 = c1(x1,x2,x3) == 0;
C2 = c2(x1,x2,x3) == 0;
C3 = c3(x1,x2,x3) == 0;
C4 = c4(x1,x2,x3) == 0;

%%
% Constraints Satisfy

x1 = 7;
x2 = 0;
x3 = 0;

C1 = c1(x1,x2,x3) <= 0
C2 = c2(x1,x2,x3) <= 0
C3 = c3(x1,x2,x3) <= 0
C4 = c4(x1,x2,x3) <= 0


% 1 Satisfy 0 Non-Satisfy

F1 = f1(x1,x2,x3) 

%%
syms x1 x2 x3 lambda1 lambda2 lambda3

d = [d_x1 
     d_x2
     d_x3];

p = [g_x1
     g_x2
     g_x3]

lagrangian = [d(1) + lambda1*p(1)
              d(2) + lambda2*p(2)
              d(3) + lambda3*p(3)]

equations = [lagrangian == 0
             constraint_active <= 0
             lambda1*(constraint_active) == 0
             lambda1 >= 0]



solution = solve(equations)

%%
d = [double(subs(d_x1,7))
     double(subs(d_x2,0))
     double(subs(d_x3,0))]

p = [double(subs(g_x1,7))
     double(subs(g_x2,0))
     double(subs(g_x3,0))]

lagrangian = sum(d)