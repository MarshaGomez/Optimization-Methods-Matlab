%% Noncooperative game theory -- Exercise 10.5

clear; close all; clc;
global C1 C2

%% data

C1 = [ 3  3
       4  1
       6  0 ];

C2 = [ 3  4
       4  0
       3  5 ] ;

[m,n] = size(C1) ;

%% check if w is a Nash equilibrium

w = [ 1/3 1/3 1/3 1/2 1/2 ]';

gap(w)
reggap(w,1)
Dgap(w,1,10)

%% find a local minimum of the regularized gap function

fprintf('Regularized gap function - local minimum\n');

alfa = 1 ;

% find a local minimum
options = optimset('Display','off');
[locmin,optval] = fmincon(@(z) reggap(z,alfa),w,[],[],...
    [ones(1,m) zeros(1,n) ; zeros(1,m) ones(1,n)], [1;1],...
    zeros(m+n,1),[],[],options)

%% try to find a global minimum of the regularized gap function
%  with a multistart approach

fprintf('Regularized gap function - multistart approach\n');

for i = 1 : 100
    
    % starting point
    x0 = rand(m+n,1);
    x0(1:m) = x0(1:m)/sum(x0(1:m));
    x0(m+1:m+n) = x0(m+1:m+n)/sum(x0(m+1:m+n));
    
    % find a local minimum
    [locmin,optval] = fmincon(@(z) reggap(z,alfa),x0,[],[],...
        [ones(1,m) zeros(1,n) ; zeros(1,m) ones(1,n)],...
        [1;1],zeros(m+n,1),ones(m+n,1),[],options);
    if optval < 1e-4
        locmin
        optval
        break
    end
end

%% try to find a global minimum of the D-gap function
%  with a multistart approach

fprintf('D-gap function - multistart approach\n');

alfa = 1 ;
beta = 10 ;

for i = 1 : 100
    
    % starting point
    x0 = rand(m+n,1);
    x0(1:m) = x0(1:m)/sum(x0(1:m));
    x0(m+1:m+n) = x0(m+1:m+n)/sum(x0(m+1:m+n));
    
    % find a local minimum
    [locmin,optval] = fminunc(@(z) Dgap(z,alfa,beta),x0,options);
    if optval < 1e-4
        locmin
        optval
        break
    end
end
