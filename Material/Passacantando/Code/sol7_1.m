%% Unconstrained optimization -- Exercise 7.1

clear; close all; clc;

%% data

Q = [6     0    -4     0
     0     6     0    -4
    -4     0     6     0
     0    -4     0     6];
c  = [ 1 -1  2 -3 ]';
x0 = [ 10  0  0  0 ]';
tolerance = 1e-6 ;

%% method

fprintf('Gradient method with exact line search\n\n');
fprintf('iter \t f(x) \t\t\t ||grad f(x)||\n\n');

iter = 0 ;

% starting point
x = x0;

while true
    v = 0.5*x'*Q*x + c'*x;
    g = Q*x + c ;
    fprintf('%2.0f \t %2.13f \t %1.9f\n',iter,v,norm(g));
    
    % stopping criterion
    if norm(g) < tolerance
        break
    end
    
    %   search direction
    d = -g;
    
    %   step size
    t = (-g'*d)/(d'*Q*d) ;
    
    %   new point
    x = x + t*d ;
    iter = iter + 1 ;
end
