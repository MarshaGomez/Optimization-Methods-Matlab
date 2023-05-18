%% Multiobjective optimization -- Exercise 9.7

close all; clear; clc;

%% data

C = [ 1  2 -3 ; 
     -1 -1 -1 ;
     -4 -2  1 ];

A = [ 1  1  1 ;
      0  0  1
      -eye(3) ];
  
b = [ 10 ; 5 ; 0 ; 0 ; 0 ] ;

%% ideal point

p = size(C,1);
n = size(C,2);
m = size(A,1);
options = optimset('Display','off');

z = zeros(p,1) ;
for i = 1 : p
    [~,z(i)] = linprog(C(i,:)',A,b,[],[],[],[],[],options);
end
z

%% goal method 

% 1-norm

gm1 = linprog([zeros(n,1);ones(p,1)], [C -eye(p); -C -eye(p); A zeros(m,p)],[z;-z;b],...
    [],[],[],[],[],options);
gm1 = gm1(1:n)


% 2-norm

gm2 = quadprog(C'*C,-C'*z,A,b,[],[],[],[],[],options)

% inf-norm

[gminf, vinf] = linprog([zeros(n,1);1], [C -ones(p,1); -C -ones(p,1); A zeros(m,1)],[z;-z;b],...
    [],[],[],[],[],options);
gminf = gminf(1:n)


% check if gminf is a minimum

c = [zeros(n,1) ; -ones(p,1)] ;
P = [C eye(p); 
     A zeros(m,p) ; 
     zeros(n,n) -eye(p)] ;
q = [C*gminf ; b ; zeros(p,1)] ;
[~,v_minimum] = linprog(c,P,q,[],[],[],[],[],options) 
