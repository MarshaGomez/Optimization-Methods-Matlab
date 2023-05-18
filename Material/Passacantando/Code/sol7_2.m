%% Unconstrained optimization -- Exercise 7.2

clear; close all; clc; 

%% data 

% the objective function is defined in f.m

alpha = 0.1;
gamma = 0.9;
tbar = 1;
x0 = [ 0 ; 0];
tolerance = 1e-3 ;

%% method

fprintf('Gradient method with Armijo inexact line search\n\n');
fprintf('iter \t f(x) \t\t ||grad f(x)||\n\n');

iter = 0 ;
x = x0 ;

while true
    [v, g] = f(x);
    fprintf('%2.0f \t %1.9f \t %1.7f\n',iter,v,norm(g));
    
    % stopping criterion
    if norm(g) < tolerance
        break
    end
    
    % search direction
    d = -g;
    
    % Armijo inexact line search
    t = tbar ;
    while f(x+t*d) > v + alpha*g'*d*t
        t = gamma*t ;
    end
        
    % new point
    x = x + t*d ;
    iter = iter + 1 ;
end
