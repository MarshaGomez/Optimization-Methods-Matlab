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

