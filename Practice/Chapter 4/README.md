# Chapter 4. Support Vector Machines for (supervised) classification problems

## 4.1 Linear SVM. Primal Model
[Exercise 4.1.](https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/Chapter%204/Exercise_4_1.mlx) Find the separating hyperplane with maximum margin for the data set given:

````matlab
close all;
clear;
clc;

matlab.lang.OnOffSwitchState = 1;

A=[ 0.4952    6.8088    
    2.6505    8.9590    
    3.4403    5.3366    
    3.4010    6.1624    
    5.2153    8.2529    
    7.6393    9.4764    
    1.5041    3.3370    
    3.9855    8.3138    
    1.8500    5.0079    
    1.2631    8.6463    
    3.8957    4.9014    
    1.9751    8.6199    
    1.2565    6.4558    
    4.3732    6.1261    
    0.4297    8.3551    
    3.6931    6.6134    
    7.8164    9.8767    
    4.8561    8.7376    
    6.7750    7.9386    
    2.3734    4.7740 
    0.8746    3.0892  
    2.3088    9.0919   
    2.5520    9.0469  
    3.3773    6.1886  
    0.8690    3.7550   
    1.8738    8.1053   
    0.9469    5.1476   
    0.9718    5.5951   
    0.4309    7.5763   
    2.2699    7.1371 ];

B=[7.2450    3.4422
   7.7030    5.0965 
   5.7670    2.8791 
   3.6610    1.5002  
   9.4633    6.5084
   9.8221    1.9383 
   8.2874    4.9380 
   5.9078    0.4489 
   4.9810    0.5962 
   5.1516    0.5319 
   8.4363    5.9467  
   8.4240    4.9696  
   7.6240    1.7988
   3.4473    0.2725 
   9.0528    4.7106  
   9.1046    3.2798   
   6.9110    0.1745  
   5.1235    3.3181  
   7.5051    3.3392  
   6.3283    4.1555   
   6.1585    1.5058
   8.3827    7.2617 
   5.2841    2.7510 
   5.1412    1.9314
   6.0290    1.9818  
   5.8863    1.0087  
   9.5110    1.3298 
   9.3170    1.0890  
   6.5170    1.4606   
   9.8621    4.3674]; 
   
% Show values distribution
xlim = ([0 10]);
ylim = ([0 10]);
plot(A(:,1), A(:,2), "r*", B(:,1), B(:,2), "b*");
````
<img src="https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/img/Chapter_4_1_1.png" width="300" height="300" />

````matlab
% Training Set
T = [A;B];

% Set Values
nA = size(A,1);
nB = size(B,1);

Q = [1 0 0
     0 1 0
     0 0 0];

D = [-A -ones(nA,1)
      B  ones(nB,1)];

d = -ones(nA+nB, 1);

f = zeros(3,1);

% Set Problem
options = optimset("LargeScale", "off", "Display", "off");
lambda = quadprog(Q,f,D,d,[],[],[],[],[],options);

% Set w and b
w = lambda(1:2);
b = lambda(3);

% Plot with the Max Margins 
x = 0:0.1:10;
u = (-w(1)/w(2)).*x - b/w(2);
v = (-w(1)/w(2)).*x + (1-b)/w(2);
vv = (-w(1)/w(2)).*x + (-1-b)/w(2);

plot(A(:,1), A(:,2), "r*", B(:,1), B(:,2), "b*");
hold on
xlim = ([0 10]);
ylim = ([0 10]);
plot(x,u, "k-", x,v,"r-", x,vv, "b-");
````

<img src="https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/img/Chapter_4_1_2.png" width="300" height="300" />

[Exercise 4.2](https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/Chapter%204/Exercise_4_2.mlx) **Linear SVM. Dual Model** Find the separating hyperplane with maximum margin for the data set given:

````matlab
A = [0.4952    6.8088    
    2.6505    8.9590    
    3.4403    5.3366    
    3.4010    6.1624    
    5.2153    8.2529    
    7.6393    9.4764    
    1.5041    3.3370    
    3.9855    8.3138    
    1.8500    5.0079    
    1.2631    8.6463    
    3.8957    4.9014    
    1.9751    8.6199    
    1.2565    6.4558    
    4.3732    6.1261    
    0.4297    8.3551    
    3.6931    6.6134    
    7.8164    9.8767    
    4.8561    8.7376    
    6.7750    7.9386    
    2.3734    4.7740 
    0.8746    3.0892  
    2.3088    9.0919   
    2.5520    9.0469  
    3.3773    6.1886  
    0.8690    3.7550   
    1.8738    8.1053   
    0.9469    5.1476   
    0.9718    5.5951   
    0.4309    7.5763   
    2.2699    7.1371];

B = [7.2450    3.4422
   7.7030    5.0965 
   5.7670    2.8791 
   3.6610    1.5002  
   9.4633    6.5084
   9.8221    1.9383 
   8.2874    4.9380 
   5.9078    0.4489 
   4.9810    0.5962 
   5.1516    0.5319 
   8.4363    5.9467  
   8.4240    4.9696  
   7.6240    1.7988
   3.4473    0.2725 
   9.0528    4.7106  
   9.1046    3.2798   
   6.9110    0.1745  
   5.1235    3.3181  
   7.5051    3.3392  
   6.3283    4.1555   
   6.1585    1.5058
   8.3827    7.2617 
   5.2841    2.7510 
   5.1412    1.9314
   6.0290    1.9818  
   5.8863    1.0087  
   9.5110    1.3298 
   9.3170    1.0890  
   6.5170    1.4606   
   9.8621    4.3674];

plot(A(:,1), A(:,2), 'r*', ...
    B(:,1), B(:,2), 'b*');
title("DataSet");

````
<img src="https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/img/Chapter_4_2_1.png" width="300" height="300" />

````matlab
% Set Environment Values
T = [A;B];
nA = size(A,1);
nB = size(B,1);

y = [ones(nA,1); -ones(nB,1)];
l = length(y);

Q = zeros(l,l);
f = -ones(l,1);
lb = zeros(l,1);

% Set Q
for i = 1:l
    for j = 1:l
        Q(i,j) = y(i)*y(j)*T(i,:)*T(j,:)';
    end 
end


% Set of Problem
options = optimset('LargeScale', 'off', 'Display', 'off');
lambda = quadprog(Q,f,[],[],y',0,lb,[],[],options);


% Get w dual
w = zeros(2,1);
for i = 1:l
    w = w + lambda(i)*y(i)*T(i,:)'; 
end

% Get b dual
indicator = find(lambda > 1e-3);
i = indicator(1);
b = 1/y(i) - w'*T(i,:)';


% Plot the solution
x = 0:0.1:10;
u = (-w(1)/w(2)).*x - b/w(2);
v = (-w(1)/w(2)).*x + (1-b)/w(2);
vv = (-w(1)/w(2)).*x + (-1-b)/w(2);

plot(A(:,1), A(:,2), 'r*', ...
     B(:,1), B(:,2), 'b*', ...
     x,u,'k-', ...
     x,v, 'r-', ...
     x,vv, 'b-');
title("Linear SVM Dual Model")

% Support Points
support = find(lambda > 1e-3);
support_A = support(support <= nA);
support_B = support(support > nA) - nA;

hold on 
plot(A(support_A,1), A(support_A,2), 'ro', ...
     B(support_B,1), B(support_B,2), 'bo', 'LineWidth',3);
xlim([0 10]);
ylim([0 10]);
hold off
````

<img src="https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/img/Chapter_4_2_2.png" width="300" height="300" />

[Exercise 4.4](https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/Chapter%204/Exercise_4_4.mlx) **Linear SVM. Dual Model. Soft Margins** Find the separating hyperplane for the data set given by solving the dual problem (5) with C = 10. What is the value of λi corresponding to the misclassified points?

````matlab
close all;
clear;
clc;
matlab.lang.OnOffSwitchState = 1;

A = [0.4952    6.8088
    2.6505    8.9590
    3.4403    5.3366
    3.4010    6.1624
    5.2153    8.2529
    7.6393    9.4764
    1.5041    3.3370
    3.9855    8.3138
    1.8500    5.0079
    1.2631    8.6463
    3.8957    4.9014
    1.9751    8.6199
    1.2565    6.4558
    4.3732    6.1261
    0.4297    8.3551
    3.6931    6.6134
    7.8164    9.8767
    4.8561    8.7376
    6.7750    7.9386
    2.3734    4.7740
    0.8746    3.0892
    2.3088    9.0919
    2.5520    9.0469
    3.3773    6.1886
    0.8690    3.7550
    1.8738    8.1053
    0.9469    5.1476
    0.9718    5.5951
    0.4309    7.5763
    2.2699    7.1371
    4.5000    2.0000];

B = [7.2450    3.4422
    7.7030    5.0965
    5.7670    2.8791
    3.6610    1.5002
    9.4633    6.5084
    9.8221    1.9383
    8.2874    4.9380
    5.9078    0.4489
    4.9810    0.5962
    5.1516    0.5319
    8.4363    5.9467
    8.4240    4.9696
    7.6240    1.7988
    3.4473    0.2725
    9.0528    4.7106
    9.1046    3.2798
    6.9110    0.1745
    5.1235    3.3181
    7.5051    3.3392
    6.3283    4.1555
    6.1585    1.5058
    8.3827    7.2617
    5.2841    2.7510
    5.1412    1.9314
    6.0290    1.9818
    5.8863    1.0087
    9.5110    1.3298
    9.3170    1.0890
    6.5170    1.4606
    9.8621    4.3674
    6.0000    8.0000
    2.0000    3.0000];

plot(A(:,1),A(:,2), 'r*',B(:,1),B(:,2), 'b*');
title('Dataset');
````

<img src="https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/img/Chapter_4_4_1.png" width="300" height="300" />

````matlab
% Set the environment values
nA = size(A,1);
nB = size(B,1);

% Training set
T = [A;B];

C = 10;
y = [ones(nA,1);-ones(nB,1)];
l = length(y);

% Set Q
Q = zeros(l,l);

for i = 1:l
    for j = 1:l
        Q(i,j) = y(i)*y(j)*T(i,:)*T(j,:)';
    end 
end 

f = -ones(l,1);
lb = zeros(l,1);
ub = C*ones(l,1);

% Set the problem
options = optimset('LargeScale', 'off', 'Display', 'off');
lambda = quadprog(Q,f,[],[],y',0, lb,ub, [], options);

% Set W
w = zeros(2,1);
for i = 1:l
    w = w + lambda(i)*y(i)*T(i,:)';
end

% Set b
indicator = find((lambda > 1e-3) & (lambda <= C-1e-3));
i = indicator(1);
b = 1/y(i) - w'*T(i,:)';

% Plot the hyperplane optimize 
x = 0:0.1:10;
u = (-w(1)/w(2)).*x - b/w(2);
v = (-w(1)/w(2)).*x + (1-b)/w(2);
vv = (-w(1)/w(2)).*x + (-1-b)/w(2);
 
plot(A(:,1), A(:,2), 'r*' ,B(:,1), B(:,2), 'b*', ...
     x, u,'k-', x,v, 'r-', x,vv, 'b-');

% Plot the support Points
support = find(lambda > 1e-3);
support_A = support(support <= nA);
support_B = support(support > nA)-nA;

hold on
plot(A(support_A,1), A(support_A,2), 'ro' ,B(support_B,1), B(support_B,2), 'bo');
xlim([0 10]);
ylim([0 10]);
hold off

````

<img src="https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/img/Chapter_4_4_2.png" width="300" height="300" />

````matlab

fprintf('Lambda Values \n Support Points A');
fprintf('%0.2f \n', lambda(support_A));

fprintf(' Support Points B');
fprintf('%0.2f \n', lambda(support_B));

````

Lambda Values


| Support Points A|  Support Points B |
|-------------------|---------------------|
| 5.01 | 0.00 |
| 10.00 | 0.00 |
| 1.40 | 10.00 |
| 10.00 | 0.00 |
| 10.00 | |

 
[Exercise 4.5.](https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/Chapter%204/Exercise_4_5.mlx) **Nonlinear SVM. Dual Model** Find the optimal separating surface for the data set given using a Gaussian kernel with parameters C = 1 and gamma = 1.

````matlab
A = [0.0113    0.2713
    0.9018   -0.1121
    0.2624   -0.2899
    0.3049    0.2100
   -0.2255   -0.7156
   -0.9497   -0.1578
   -0.6318    0.4516
   -0.2593    0.6831
    0.4685    0.1421
   -0.4694    0.8492
   -0.5525   -0.2529
   -0.8250    0.2802
    0.4463   -0.3051
    0.3212   -0.2323
    0.2547   -0.9567
    0.4917    0.6262
   -0.2334    0.2346
    0.1510    0.0601
   -0.4499   -0.5027
   -0.0967   -0.5446];

B = [1.2178    1.9444
  -1.8800    0.1427
  -1.6517    1.2084
   1.9566   -1.7322
   1.7576   -1.9273
   0.7354    1.1349
   0.1366    1.5414
   1.5960    0.5038
  -1.4485   -1.1288
  -1.2714   -1.8327
  -1.5722    0.4658
   1.7586   -0.5822
  -0.3575    1.9374
   1.7823    0.7066
   1.9532    1.0673
  -1.0233   -0.8180
   1.0021    0.3341
   0.0473   -1.6696
   0.8783    1.9846
  -0.5819    1.8850];

% Plot the Dataset
plot(A(:,1), A(:,2), 'r*', B(:,1), B(:,2), 'b*');
````

<img src="https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/img/Chapter_4_5_1.png" width="300" height="300" />

````matlab
% Set environment variables
gamma = 1;
C = 1;
nA = size(A,1);
nB = size(B,1);

% Training Set 
T = [A;B];

y = [ones(nA,1); -ones(nB,1)];
l = length(y);

% Quadratic variables
f = -ones(l,1);
lb = zeros(l,1);
ub = C*ones(l,1);

K = zeros(l,l);
for i = 1:l
    for j = 1:l
        K(i,j) = exp(-gamma*norm(T(i,:)-T(j,:))^2);
    end
end

% Compute Q
Q = zeros(l,l);
for i = 1:l
    for j = 1:l
        Q(i,j) = y(i)*y(j)*K(i,j);
    end 
end

% Set the Problem
options = optimset('Largescale', 'off', 'Display', 'off');
[lambda, value] = quadprog(Q,f,[],[],y',0,lb,ub,[],options);

% Get b
indicator = find((lambda > 1e-3) & (lambda < C-1e-3));
i = indicator(1);
b = 1/y(i);
for j = 1:l
    b = b - lambda(j)*y(j)*K(i,j);
end

% Get general points 
AA = [];
BB = [];

for xx = -2:0.1:2
    for yy = -2:0.1:2
        s = 0;
        for i = 1:l 
            s = s + lambda(i)*y(i)*exp(-gamma*norm(T(i,:) - [xx yy])^2);
        end
        s = s + b;
        if s > 0 
            AA = [AA; xx yy];
        else 
            BB = [BB; xx yy];
        end
    end
end

% Plot the Nonlinear Support Vector Machine Dual Model Graph
plot(A(:,1), A(:,2), 'ro', B(:,1), B(:,2), 'bo', 'LineWidth',5);
hold on
plot(AA(:,1), AA(:,2), 'r.', BB(:,1), BB(:,2), 'b.');
hold off
````

<img src="https://github.com/MarshaGomez/Optimization-Matlab-Exams/blob/master/Practice/img/Chapter_4_5_2.png" width="300" height="300" />

