# Chapter 9

# Multiobjective optimization

[Exercise 9.1.](https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/Chapter%209/Exercise_9_1.mlx) **Linear Programming. Multiobjective optimization** Consider the linear multiobjective problem:

<img src="https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/img/Chapter_9_1_2.png" />

Check if the points u = (5, 0, 5), v = (4, 4, 2) and w = (1, 4, 4) are minima or weak minima by solving the corresponding auxiliary problems.

````matlab
close all;
clear;
clc;

matlab.lang.OnOffSwitchState = 1;

x1 = -10:0.1:10;
x2 = -10:0.1:10;

[X1,X2] = meshgrid(x1,x2);

X3 = (X1+2.*X2)/3;
X4 = (-X1-X2);
X5 = (4.*X1+2.*X2);
surface(X1,X2,X3, 'FaceAlpha',0.5);
hold on
surface(X1,X2,X4,'FaceAlpha',0.2);
surface(X1,X2,X5, 'FaceAlpha',0.7);
view(3);
hold off
````

<img src="https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/img/Chapter_9_1_1.png" width="300" height="300" />

````matlab
% Environment Variables
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

u = [5
     0
     5];

v = [4
     4
     2];

w = [1
     4
     4];

y = v;

% Variables number
n = size(C,2);
% Functions number
p = size(C,1);
% Constraints number
m = size(A,1);

% Calculate the linear Problem

c = [zeros(n,1)
    -ones(p,1)];

P = [C           eye(p)
     A           zeros(m,n)
     zeros(p,n) -eye(p)];

q = [C*y
     b
     zeros(p,1)];

% Solve the Linear Programming check minimum
options = optimset('Display', 'off');
[~, minimum] = linprog(c,P,q, [],[],[],[],options)

% Solve the Linear Programming check weak minimum

c = [zeros(n,1)
     zeros(p,1)
     -1];

P = [zeros(p,n)  -eye(p)    ones(p,1)
     C           eye(p)     zeros(p,1)
     A           zeros(m,n) zeros(m,1)
     zeros(p,n) -eye(p)     zeros(p,1)];

q = [zeros(p,1)
     C*y
     b
     zeros(p,1)];

[~, minimum_weak] = linprog(c,P,q, [],[],[],[],options)
````

| Point | Minimum | Weak Minimum |
|-------|---------|--------------|
| (4,4,2) | -13.00 | 0.00  |
| (5,0,5) | 0.00   | 0.00  |
| (1,4,4) | -16.75 | -1.00 |


[Exercise 9.4.](https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/Chapter%209/Exercise_9_4.mlx) **Scalarization Method. Multiobjective optimization** Consider the linear multiobjective problem:

<img src="https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/img/Chapter_9_4_6.png" />

Find the set of minima and weak minima by means of the scalarization method:

````matlab
close all;
clear;
clc;

matlab.lang.OnOffSwitchState = 1;


% Plot Graph

x1 = -10:0.1:10;
x2 = -10:0.1:10;

[X1,X2] = meshgrid(x1,x2);

X3 = X1-X2;
X4 = X1+X2;

surface(X1,X2,X3);
hold on
surface(X1,X2,X4);
view(3)
hold off
````


<img src="https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/img/Chapter_9_4_1.png" width="300" height="300" />

````matlab
% 2D plot X1
x1 = x2;
plot(x1,x2);
````

<img src="https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/img/Chapter_9_4_2.png" width="300" height="300" />

````matlab
% 2D plot X2
x2 = -x1;
plot(x1,x2);
````

<img src="https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/img/Chapter_9_4_3.png" width="300" height="300" />

````matlab
% Feasible Area
x1 = -10:0.1:10;
x2 = -10:0.1:10;

[X1, X2] = meshgrid(x1,x2);

condition_x1 = X1((X2 <= 2*X1) & (-X1 - X2 <= 0) & (5*X1 - X2 <= 6));
condition_x2 = X2((X2 <= 2*X1) & (-X1 - X2 <= 0) & (5*X1 - X2 <= 6));

plot(condition_x1, condition_x2, '.');
hold on
plot([0 2],[0 4],'-k');
plot([2 1],[4 -1],'-k');
plot([1 0],[-1 0],'-k');

````

<img src="https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/img/Chapter_9_4_4.png" width="300" height="300" />


````matlab
% Linear Multiobjective problem. Scalarization Method

C = [1 -1
     1  1];

A = [-2  1
     -1 -1
      5 -1];

b = [0
     0
     6];

% Calculate different Alpha
options = optimset('Display', 'off');
for alpha = 0.001:0.001:0.999
    f = [1
         1-2*alpha];
  
    x = linprog(f,A,b, [], [], [],[], options);
    plot(x(1),x(2), 'r.');
end

fprintf('Alpha 0 \n');
alpha = 0;
f = [1
     1-2*alpha];

x0 = linprog(f,A,b, [], [], [],[], options);
plot(x0(1), x0(2), 'ko');


alpha = 1;

f = [1
    1-2*alpha];

x1 = linprog(f,A,b, [], [], [], [], options);
plot(x1(1), x1(2), 'bo');
hold off
````

<img src="https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/img/Chapter_9_4_5.png" width="300" height="300" />


[Exercise 9.5.](https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/Chapter%209/Exercise_9_5.mlx) **Scalarization Method. Multiobjective optimization. Quadratic Problem. Nonlinear Problem** Consider the nonlinear multiobjective problem:

<img src="https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/img/Chapter_9_5_1.png"/>

* a) Find the set of weak minima by means of the scalarization method. 
* b) What is the set of minima?


````matlab
close all;
clear;
clc;

matlab.lang.OnOffSwitchState = 1;

x1 = 0:0.1:10;
x2 = 0:0.1:10;

[X1, X2] = meshgrid(x1,x2);

X3 = X1;
X4 = X1.^2 + X2.^2 -2.*X1;

surface(X1,X2,X3);
hold on
surface(X1,X2,X4);
view(3);
hold off
````

<img src="https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/img/Chapter_9_5_2.png" width="300" height="300" />


````matlab
% 2D
x4 = sqrt(-x1.^2 + 2.*x1);
x3 = x1;

plot(x1,x4);
hold on
plot(x1,x3);
hold off
````

<img src="https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/img/Chapter_9_5_3.png" width="300" height="300" />



````matlab
condition_x1 = X1((-X1<=0) & (-X2<=0) & (X1+X2<=2)); 
condition_x2 = X2((-X1<=0) & (-X2<=0) & (X1+X2<=2));

plot(condition_x1, condition_x2, 'b.');

hold on

Q = [2 0
     0 2];

c = [-2
      0];

A = [ -1  0
       0 -1
       1  1];

b = [0
     0
     2];

options = optimset('Display', 'off');

for alpha = 0.001:0.001:0.999
    Q_alpha = Q*alpha;
    c_alpha = [1-3*alpha ; 0];
    x = quadprog(Q_alpha, c_alpha, A,b,[],[],[],[],[],options);
    plot(x(1), x(2), 'r.');
end

% Alpha 0
alpha = 0;
Q_alpha = Q*alpha;
c_alpha = [1-3*alpha ; 0];
x = quadprog(Q_alpha, c_alpha, A,b,[],[],[],[],[],options);
plot(x(1), x(2), 'ro');

% Alpha 1
alpha = 1;
Q_alpha = Q*alpha;
c_alpha = [1-3*alpha ; 0];
x = quadprog(Q_alpha, c_alpha, A,b,[],[],[],[],[],options);
plot(x(1), x(2), 'go');

hold off
````

<img src="https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/img/Chapter_9_5_4.png" width="300" height="300" />

[Exercise 9.6.](https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/Chapter%209/Exercise_9_6.mlx) **Noninear Programming. Multiobjective optimization. Scalarization Method** Consider the nonlinear multiobjective problem:

<img src="https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/img/Chapter_9_6_4.png" />

Find the set of minima and weak minima by means of the scalarization method.

````matlab
close all;
clear;
clc;

matlab.lang.OnOffSwitchState = 1;

% Plot of data 
x1 = -10:0.1:10;
x2 = -10:0.1:10;


[X1,X2] = meshgrid(x1,x2);

X3 = X1.^2 + X2.^2 + 2.*X1 - 4.*X2;
X4 = X1.^2 + X2.^2 - 6.*X1 - 4.*X2;

surface(X1,X2,X3);
hold on
surface(X1,X2,X4);
view(3);
hold off
````

<img src="https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/img/Chapter_9_6_1.png" width="300" height="300"  />


````matlab
% Check convexity
Q1 = [2 0
      0 2];

f1 = [2
     -4];

Q2 = [2 0
      0 2];

f2 = [-6
      -4];

eigenvalue_1 = eig(Q1)
eigenvalue_2 = eig(Q1)
````

eigenvalue_1 | eigenvalue_2
----|----
2 | 2 
2 | 2


````matlab
% Feasible region
x1 = -10:0.02:10;
x2 = -10:0.02:10;


[X1,X2] = meshgrid(x1,x2);

condition_x1 = X1((-X2<=0) & (-2.*X1 + X2<=0) & (2.*X1+X2<=4));
condition_x2 = X2((-X2<=0) & (-2.*X1 + X2<=0) & (2.*X1+X2<=4));

plot(condition_x1, condition_x2, 'b.');
````

<img src="https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/img/Chapter_9_6_2.png" width="300" height="300" />

````matlab
% Set the environment values
A = [0 -1
    -2  1
     2  1];

b = [0
     0
     4];

% Set the problem
options = optimset('Display', 'off');

hold on
for alpha = 0.001:0.001:0.999
    f = [8*alpha-6
         -4];
    x = quadprog(Q1,f, A, b, [], [], [], [], [], options);
    plot(x(1), x(2), 'r.');
end

% Extreme alpha
% alpha = 0;
x0 = quadprog(Q2, f2, A, b, [], [], [], [], [], options);
plot(x0(1), x0(2), 'go');

%alpha = 1;
x1 = quadprog(Q1, f1, A, b, [], [], [], [], [], options);
plot(x1(1), x1(2), 'ko');

hold off
````

<img src="https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/img/Chapter_9_6_3.png" width="300" height="300" />




