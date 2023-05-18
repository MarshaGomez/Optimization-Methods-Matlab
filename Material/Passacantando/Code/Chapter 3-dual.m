clear; close all; clc;

%% Primal problem

x = -1 : 0.01 : 3;
y = 2*x.^4+x.^3-20*x.^2+x;
[xopt,vp]=fminbnd(@(x) 2*x^4+x^3-20*x^2+x,-1,3);
plot(x,y,'b-',x,vp*ones(length(x),1),'b-','LineWidth',1.5);
title('Primal problem:      min 2x^4+x^3-20x^2+x    s.t.    x^2-2x-3 \leq 0');
axis([-1 3 -45 15]);

pause

%% Dual problem

x = -4 : 0.01 : 4;
phi = [];
figure;
for lam = 0 : 0.01 : 5
    y = 2*x.^4+x.^3-20*x.^2+x + lam*(x.^2-2.*x-3);
    
    [~,vl1]=fminbnd(@(x) 2*x^4+x^3-20*x^2+x+lam*(x^2-2*x-3),-3,0);
    [~,vl2]=fminbnd(@(x) 2*x^4+x^3-20*x^2+x+lam*(x^2-2*x-3),0,3);
    vl = min(vl1,vl2);
    phi = [phi ; vl];
    plot(x,y,'r-',x,vl*ones(length(x),1),'r-',x,vp*ones(length(x),1),'b-','LineWidth',1.5);
    title(['Lagrangian function with \lambda = ',num2str(lam)]);
    axis([-4 4 -70 50]);
    if lam == 0
        pause
    else
        pause(0.03)
    end
end

vd = max(phi);
figure;
plot(0:0.01:5,phi,'r-',0:0.01:5,vd*ones(length(0:0.01:5),1),'r-',...
    0:0.01:5,vp*ones(length(0:0.01:5),1),'b-','LineWidth',1.5);
title('Dual problem');
