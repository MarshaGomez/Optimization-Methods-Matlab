function v = reggap(z,alfa)

% regularized gap function of a bimatrix game

global C1 C2;

[m,n] = size(C1) ;
x = z(1:m) ;
y = z(m+1:m+n) ;

options = optimset('Display','off');

[~,v1] = quadprog(alfa*eye(m),C1*y-alfa*x,[],[],ones(1,m),1,...
    zeros(m,1),ones(m,1),[],options);

[~,v2] = quadprog(alfa*eye(n),C2'*x-alfa*y,[],[],ones(1,n),1,...
    zeros(n,1),ones(n,1),[],options);

v = x'*(C1+C2)*y - 0.5*alfa*(norm(x)^2 + norm(y)^2) - v1 - v2 ;

end