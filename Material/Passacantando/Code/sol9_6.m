%% Multiobjective optimization -- Exercise 9.6

close all; clear; clc; hold on

%% data

A = [ 0 -1 ; 
     -2  1 ;
      2  1 ];
b = [ 0 ; 0 ; 4 ] ;

%% plot the feasible region

plot([0 2],[0 0],'-k');
plot([2 1],[0 2],'-k');
plot([1 0],[2 0],'-k');

%% solve the scalarized problem with 0 <= alfa <= 1

options = optimset('Display','off');
for alfa = 0 : 0.001 : 1
    x = quadprog([2 0 ; 0 2],[8*alfa-6 ; -4],A,b,...
        [],[],[],[],[],options) ;
    plot(x(1),x(2),'g.');
end
