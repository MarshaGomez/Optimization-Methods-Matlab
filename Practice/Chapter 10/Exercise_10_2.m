%% 
% *Game Theory*
% 
% 
% 
% 
% 
% Are there strictly dominated strategies? Are there pure strategies Nash equilibria? 
% Find a mixed strategies Nash equilibrium.

close all;
clear;
clc;

format rational;

% Game Theory

% Cost Matrix

C = [ 1  2  3
      3 -2  3
      1  2  1];

% Get dimension of Cost
[m,n] = size(C);

f = [zeros(m,1)
     1];

A = [C' -ones(n,1)];
b = [zeros(n,1)];

Aeq = [ones(1,m) 0];
beq = 1;
lb = [zeros(m,1) 
     -inf];
up = [];

options = optimset('Display', 'off');
[solution,v, ~, ~, lambda] = linprog(f,A,b,Aeq,beq,lb,up,options);

% Nash equilibrium

x = solution(1:m)
y = lambda.ineqlin