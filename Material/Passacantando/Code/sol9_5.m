%% Multiobjective optimization -- Exercise 9.5

close all; clear; clc; hold on

%% data

A = [ -1   0 ; 
       0  -1 ;
       1   1 ];
b = [ 0 ; 0 ; 2 ] ;

%% plot the feasible region

plot([0 0],[0 2],'-k');
plot([0 2],[2 0],'-k');
plot([2 0],[0 0],'-k');

%% solve the scalarized problem with 0 < alfa < 1

options = optimset('Display','off');
for alfa = 0.001 : 0.001 : 0.999
    x = quadprog([2*alfa 0 ; 0 2*alfa],[1-3*alfa ; 0],A,b,...
        [],[],[],[],[],options) ;
    plot(x(1),x(2),'g.');
end

%% solve the scalarized problem with alfa = 0

alfa = 0;
x0 = quadprog([2*alfa 0 ; 0 2*alfa],[1-3*alfa ; 0],A,b,...
    [],[],[],[],[],options) ;
plot(x0(1),x0(2),'ro');

%% solve the scalarized problem with alfa = 1

alfa = 1;
x1 = quadprog([2*alfa 0 ; 0 2*alfa],[1-3*alfa ; 0],A,b,...
        [],[],[],[],[],options) ;
plot(x1(1),x1(2),'bo');

