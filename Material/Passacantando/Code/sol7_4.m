%% Unconstrained optimization -- Exercise 7.4

clear; close all; clc;

%% data

% the objective function is defined in f.m

alpha = 0.1;
gamma = 0.9;
tbar = 1;
x0 = [ 0 0 ]';
tolerance = 1e-3 ;

%% method

fprintf('Newton method with line search\n\n');
fprintf('iter \t f(x) \t\t ||grad f(x)||\n\n');

iter = 0 ;
x = x0 ;

while true
    [v, g, H] = f(x);
    fprintf('%1.0f \t %1.7f \t %1.4e\n',iter,v,norm(g));
    
    % stopping criterion
    if norm(g) < tolerance
        break
    end
    
    % search direction H*d = -g
    d = -H\g;
    
    % Armijo inexact line search
    t = tbar ;
    while f(x+t*d) > v + alpha*g'*d*t
        t = gamma*t ;
    end
        
    % new point
    x = x + t*d ;
    iter = iter + 1 ;
    
end

