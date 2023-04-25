%% 
% *Merit Function. Game Theory*
% 
% Consider the following bimatrix game:
% 
% $$C_1 =\left\lbrack \begin{array}{cc}3 & 3\\4 & 1\\6 & 0\end{array}\right\rbrack 
% \;\;\;C_2 =\left\lbrack \begin{array}{cc}3 & 4\\4 & 0\\3 & 5\end{array}\right\rbrack$$
% 
% a) Implement in MATLAB the gap function, the regularized gap function and 
% the D-gap function. 
% 
% b) Exploit the gap function ψ to check if the point w = (x, y), where x = 
% (1/3, 1/3, 1/3) and y = (1/2, 1/2), is a Nash equilibrium. 
% 
% c) Find a local minimum of the regularized gap function ψα with α = 1 starting 
% from w. 
% 
% d) Try to find a global minimum of the regularized gap function ψα with a 
% multistart approach. 
% 
% e) Try to find a global minimum of the D-gap function ψα,β, with α = 1 and 
% β = 10, with a multistart approach.

% Cost Function. Bimatrix Game

global C1 C2

C1 = [3 3
      4 1
      6 0];

C2 = [3 4
      4 0
      3 5];

[m,n] = size(C1);

%% check if w is a Nash equilibrium

w = [ 1/3 1/3 1/3 1/2 1/2 ]';

v1 = gap(w);
v2 = reggap(w,1);
v3 = dgap(w,1,10);

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
    [locmin,optval] = fminunc(@(z) dgap(z,alfa,beta),x0,options);
    if optval < 1e-4
        locmin
        optval
        break
    end
end
%% 
% *Implementation Functions*

function v = gap(z)
    % GAP function of bimatrix game
    global C1 C2

    [m,n] = size(C1);

    x = z(1:m); 
    y = z(m+1:m+n);

    % Solution with linear programming
    options = optimset('Display', 'off');

    [~,v1] = linprog(C1*y, [],[], ones(1,m), 1, zeros(m,1), [], options);
    [~,v2] = linprog(C2'*x, [], [], ones(1,n), 1, zeros(n,1), [], options);

    v = x'*(C1+C2)*y - v1 - v2;

    % Solution without linear programming

    % v = x'*(C1+C2)*y - min(C1*y) - min(C2'*x); 
end

function v = dgap(z,alpha,beta)
    % DGAP function of a bimatrix game
    v = reggap(z, alpha) - reggap(z,beta); 
end

function v = reggap(z,alpha)

    global C1 C2
    % Regularized GAP function of a bimatrix game 

    [m,n] = size(C1);
    x = z(1:m);
    y = z(m+1:m+n);

    % Solution with Quadratic programming
    H = [alpha*eye(m)];
    f = [C1*y-alpha*x];

    Aeq = ones(1,m);
    beq = 1;
    lb = zeros(m,1);
    ub = ones(m,1);

    options = optimset('Display', 'off');
    [~, v1] = quadprog(H,f,[],[],Aeq,beq,lb,ub,[],options);
    
    H = [alpha*eye(n)];
    f = [C2'*x-alpha*y];

    Aeq = ones(1,n);
    beq = 1;
    lb = zeros(n,1);
    ub = ones(n,1);
    [~, v2] = quadprog(H,f,[],[],Aeq,beq,lb,ub,[],options);

    v = x'*(C1+C2)*y - (1/2)*alpha*(norm(x)^2 + norm(y)^2) - v1 - v2;
end