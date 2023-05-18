    %% Multiobjective optimization -- Exercise 9.4

close all; clear; clc; hold on

%% data

A = [ -2   1 ; 
      -1  -1 ;
       5  -1 ];
b = [ 0 ; 0 ; 6 ] ;

%% plot the feasible region

plot([0 2],[0 4],'-k');
plot([2 1],[4 -1],'-k');
plot([1 0],[-1 0],'-k');

%% solve the scalarized problem with 0 < alfa < 1

options = optimset('Display','off');
for alfa = 0.001 : 0.001 : 0.999
    x = linprog([1 ; 1-2*alfa],A,b,[],[],[],[],options) ;
    plot(x(1),x(2),'g.');
end

%% solve the scalarized problem with alfa = 0

alfa = 0;
x0 = linprog([1 ; 1-2*alfa],A,b,[],[],[],[],options) ;
plot(x0(1),x0(2),'ro');

%% solve the scalarized problem with alfa = 1

alfa = 1;
x1 = linprog([1 ; 1-2*alfa],A,b,[],[],[],[],options) ;
plot(x1(1),x1(2),'bo');

