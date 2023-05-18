%% Unconstrained optimization -- Exercise 7.3

clear; close all; clc;

%% data 

Q = [6     0    -4     0
     0     6     0    -4
    -4     0     6     0
     0    -4     0     6];
c = [ 1 -1  2 -3 ]';
x0 = [ 0 0 0 0 ]';
tolerance = 1e-6 ;

%% method

fprintf('Conjugate Gradient method\n\n');
fprintf('iter \t f(x) \t\t ||grad f(x)||\n\n');

iter = 0;

% starting point
x = x0;

while true
    v = 0.5*x'*Q*x + c'*x;
    g = Q*x + c ;
    fprintf('%1.0f \t %1.4f \t %1.4e\n',iter,v,norm(g));
    
    % stopping criterion
    if norm(g) < tolerance
        break
    end
    
    %   search direction
    if iter == 0
        d = -g; 
    else
        beta = (norm(g)^2)/(norm(g_prev)^2);
        d = -g + beta*d_prev;
    end
    
    %   step size
    t = (norm(g)^2)/(d'*Q*d);
    
    %   new point
    iter = iter + 1;
    x = x + t*d;
    d_prev = d ;
    g_prev = g ; 
end

