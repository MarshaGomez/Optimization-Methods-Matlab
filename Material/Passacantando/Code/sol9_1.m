%% Multiobjective optimization -- Exercise 9.1

close all; clear; clc;

%% data

C = [ 1  2 -3 ;
     -1 -1 -1 ;
     -4 -2  1 ];

A = [ 1  1  1 ;
      0  0  1 ;
      -eye(3) ];

b = [ 10 ; 5 ; 0 ; 0 ; 0 ] ;

% given point

% y = [ 5 ; 0 ; 5 ] ;
% y = [ 4 ; 4 ; 2 ] ;
y = [ 1 ; 4 ; 4 ] ;

%% solve the problem

n = size(C,2);
p = size(C,1);
m = size(A,1);

% check if y is a minimum

c = [zeros(n,1) ; -ones(p,1)] ;
P = [C eye(p);
    A zeros(m,p) ;
    zeros(n,n) -eye(p)] ;
q = [C*y ; b ; zeros(p,1)] ;
options = optimset('Display','off');
[~,v_minimum] = linprog(c,P,q,[],[],[],[],[],options)

% check if y is a weak minimum

c = [ zeros(n,1) ; zeros(p,1) ; -1 ] ;
P = [zeros(p,n) -eye(p) ones(p,1) ;
    C eye(p) zeros(p,1) ;
    A zeros(m,p) zeros(m,1);
    zeros(n,n) -eye(p) zeros(p,1) ] ;
q = [zeros(p,1) ; C*y ; b ; zeros(p,1)] ;
[~,v_weak_minimum] = linprog(c,P,q,[],[],[],[],[],options)
