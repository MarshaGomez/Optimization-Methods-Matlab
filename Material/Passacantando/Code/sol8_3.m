%% Constrained optimization -- Exercise 8.3

clear; close all; clc;

%% data 

global Q c A b eps;

Q = [ 1 0 ; 0 2 ] ;
c = [ -3 ; -4 ] ;
A = [-2 1 ; 1 1 ; 0 -1 ];
b = [ 0 ; 4 ; 0 ];

delta = 1e-6 ;
tau = 0.1 ;
eps1 = 1 ;
x0 = [ 1 ; 1 ];

%% method

fprintf('Logarithmic barrier method\n\n');
fprintf('eps \t\t x(1) \t\t x(2) \t\t gap \n\n');

options = optimoptions('fminunc','GradObj','on',...
    'Algorithm','quasi-newton','Display','off');

x = x0;
eps = eps1 ;
m = size(A,1) ;

while true
    x = fminunc(@logbar,x,options);
    gap = m*eps;
    fprintf('%1.2e \t %1.6f \t %1.6f \t %1.2e\n',eps,x(1),x(2),gap);
    if gap < delta
        break
    else
        eps = eps*tau;
    end
end

%% logarithmic barrier function

function [v,g] = logbar(x)

    global Q c A b eps

    v = 0.5*x'*Q*x + c'*x ;
    g = Q*x + c ;

    for i = 1 : length(b)
        v = v - eps*log(b(i)-A(i,:)*x) ;
        g = g + (eps/(b(i)-A(i,:)*x))*A(i,:)' ;
    end

end