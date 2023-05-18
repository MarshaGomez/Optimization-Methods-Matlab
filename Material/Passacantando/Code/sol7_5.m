%% Unconstrained optimization -- Exercise 7.5

clear; close all; clc;

%% data

% the objective function is defined in f.m

x0 = [ 0 ; 0 ];
t0 = 5 ;
beta = 0.5 ;
epsilon = 1e-5 ;

D = [ 1 0 -1 0 ; 
      0 1 0 -1 ] ;  
  
% D = [ 1 0 -1 ; 
%       0 1 -1 ] ;

%% method

fprintf('Directional direct-search method\n\n');
fprintf(' x(1) \t\t x(2) \t\t f(x) \n\n');

x = x0;
t = t0;

v = f(x) ;
fprintf('%1.6f \t %1.6f \t %1.6f\n',x(1),x(2),v);
plot(x(1), x(2),'r.');
axis([-6 6 -6 6])
hold on ;
pause
iter = 0;

while t > epsilon   
    iter = iter + 1 ;
    newv = v ;
    i = 0 ;
    while (newv >= v) && (i < size(D,2))
        i = i + 1 ;
        newx = x+t*D(:,i);
        newv = f(newx) ;
        if newv >= v 
            plot(newx(1),newx(2),'bs');
            pause
        end
    end
    if newv < v
        x = newx ;
        v = newv ;
        plot(x(1), x(2),'r.');
        fprintf('%1.6f \t %1.6f \t %1.6f\n',x(1),x(2),v);
        pause
    else
        t = beta*t ;
    end   
end

