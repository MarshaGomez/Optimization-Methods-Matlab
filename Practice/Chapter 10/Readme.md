# Chapter 10. Noncooperative Game Theory 

## Penalty Kick
[Exercise 10.2.](https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/Chapter%2010/Exercise_10_2.mlx) Consider the following matrix game::

<table>
  <tr>
    <th> </th>
    <th colspan="4">Player 2</th>
  </tr>
  <tr>
    <th> </th>
    <th> </th>
    <th>L</th>
    <th>C</th>
    <th>R</th>
  </tr>
  <tr>
    <th rowspan="4">Player 1</th>
  </tr>
   <tr>
    <th>L</th>
    <td>1</td>
    <td>2</td>
    <td>3</td>
  </tr>
   <tr>
    <th>C</th>
    <td>3</td>
    <td>-2</td>
    <td>3</td>
  </tr>
   <tr>
    <th>R</th>
    <td>3</td>
    <td>2</td>
    <td>1</td>
  </tr>
</table>



````matlab
% Cost Matrix

C = [ 1  2  3
      3 -2  3
      1  2  1];

% Get dimension of Cost
[m,n] = size(C);

f = [zeros(m,1)
     1];

A = [C' -ones(n,1)];
b = [zeros(n,1)];

Aeq = [ones(1,m) 0];
beq = 1;
lb = [zeros(m,1) 
     -inf];
up = [];

options = optimset('Display', 'off');
[solution,v, ~, ~, lambda] = linprog(f,A,b,Aeq,beq,lb,up,options);

% Nash equilibrium

x = solution(1:m)
y = lambda.ineqlin
````

## Merit Function
[Exercise 10.5.](https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/Chapter%2010/Exercise_10_5.mlx) Consider the following bimatrix game:


$$ C_{1} = \begin{pmatrix}
3 & 3\\ 
4 & 1\\ 
6 & 0
\end{pmatrix}   C_{2} = \begin{pmatrix}
3 & 4\\ 
4 & 0\\ 
3 & 5
\end{pmatrix} $$

* a) Implement in MATLAB the gap function, the regularized gap function and the D-gap function.
* b) Exploit the gap function ψ to check if the point w = (x, y), where x = (1/3, 1/3, 1/3) and y = (1/2, 1/2), is a Nash equilibrium.
* c) Find a local minimum of the regularized gap function ψα with α = 1 starting from w.
* d) Try to find a global minimum of the regularized gap function ψα with a multistart approach.
* e) Try to find a global minimum of the D-gap function ψα,β, with α = 1 and β = 10, with a multistart approach.


````matlab
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
````
