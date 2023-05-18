%% Constrained optimization -- Exercise 8.2

clear; close all; clc;

%% data 

global Q c A b eps;

%% data

Q = [ 1 0 ; 0 2 ] ;
c = [ -3 ; -4 ] ;
A = [-2 1 ; 1 1 ; 0 -1 ];
b = [ 0 ; 4 ; 0 ];

tau = 0.1 ;
eps0 = 5 ;
tolerance = 1e-6 ;

%% method

fprintf('Penalty method\n\n');
fprintf('iter \t eps \t\t x(1) \t\t x(2) \t\t max(Ax-b)\n\n');

options = optimoptions('fminunc','GradObj','on',...
    'Algorithm','quasi-newton','Display','off');

eps = eps0;
x = [0;0];
iter = 0;

while true
    x = fminunc(@p_eps,x,options);      
    infeas = max(A*x-b);
    fprintf('%2.0f \t %1.2e \t %1.6f \t %1.6f \t %1.3e\n',iter,eps,x(1),x(2),infeas);
    if infeas < tolerance
        break
    else
        eps = tau*eps;
        iter = iter + 1 ;
    end
end

%% penalized function 

function [v,g] = p_eps(x) 

    global Q c A b eps;

    v = 0.5*x'*Q*x + c'*x ;
    g = Q*x + c ;

    for i = 1 : size(A,1)
        v = v + (1/eps)*(max(0,A(i,:)*x-b(i)))^2 ;
        g = g + (2/eps)*(max(0,A(i,:)*x-b(i)))*A(i,:)';
    end

end
