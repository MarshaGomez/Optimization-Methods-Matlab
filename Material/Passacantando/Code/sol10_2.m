%% Noncooperative game theory -- Exercise 10.2

clear; close all; clc; 
format rat;

%% cost matrix

C = [ 1  2 3
      3 -1 3
      3  2 1];

%% solve the LP problem

[m,n] = size(C) ;

[sol,v,~,~,lambda] = linprog([zeros(m,1);1],...
    [C' -ones(n,1)], zeros(n,1),...
    [ones(1,m) 0], 1,...
    [zeros(m,1); -inf], []);

% Nash equilibrium

x = sol(1:m)
y = lambda.ineqlin


