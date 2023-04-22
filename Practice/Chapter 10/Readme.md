# Chapter 10

## Noncooperative Game Theory 

[Exercise 10.2.](https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/Chapter%2010/Exercise_10_2.mlx) **Penalty kick** Consider the following matrix game::

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
