function [ p , varargout ] =  LDBCQP( BCQP , varargin )

%function [ p , status , l , v , x ] =
%                              LDBCQP( BCQP , dolh , eps , MaxFeval , m1 ,
%                                      m2 , astart , tau , sfgrd , mina )
%
% Solve the convex Box-Constrained Quadratic program
%
%  (P) min { (1/2) x^T * Q * x + q * x : 0 <= x <= u }
%
% encoded in the structure BCQP using a dual approach, where Q must be
% strictly positive definite. The box constraints 0 <= x <= u are relaxed
% (with Lagrangian multipliers \lambda^- and \lambda^+, respectively) and
% the corresponding Lagrangian Dual is solved by means of an ad-hoc
% implementation of the Projected Gradient method (since the Lagrangian
% multipliers are constrained in sign, and the Lagrangian function is
% differentiable owing to strict-positive-definiteness of Q) using a
% classical Armijo-Wolfe line search.
%
% A rough Lagrangian heuristic is implemented whereby the dual solution is
% projected on the box at each iteration to provide an upper bound, which
% is then used in the stopping criterion.
%
% Input:
%
% - BCQP, the structure encoding the BCQP to be solved within its fields:
%
%   = BCQP.Q: n \times n symmetric positive semidefinite real matrix
%
%   = BCQP.q: n \times 1 real vector
%
%   = BCQP.u: n \times 1 real vector > 0
%
% - dolh (logical scalar, optional, default value 1): true if a quick
%   Lagrangian Heuristic is implemented whereby the current solution of the
%   Lagrangian subproblem is projected on the box constraints
%
% - eps (real scalar, optional, default value 1e-6): the accuracy in the 
%   stopping criterion. This depends on the dolh parameter, see above: if
%   dolh = true then the algorithm is stopped when the relative gap
%   between the value of the best dual solution (the current one) and the
%   value of the best upper bound obtained so far (by means of the
%   heuristic) is less than or equal to eps, otherwise the stopping
%   criterion is on the (relative, w.r.t. the first one) norm of the
%   (projected) gradient
%
% - MaxFeval (integer scalar, optional, default value 1000): the maximum
%   number of function evaluations (hence, iterations will be not more than
%   MaxFeval because at each iteration at least a function evaluation is
%   performed, possibly more due to the line search).
%
% - m1 (real scalar, optional, default value 0.01): first parameter of the
%   Armijo-Wolfe-type line search (sufficient decrease). Has to be in (0,1)
%
% - m2 (real scalar, optional, default value 0.9): the second parameter of
%   the Armijo-Wolfe-type line search (strong curvature condition), it
%   should to be in (0,1);
%
% - astart (real scalar, optional, default value 1): starting value of
%   alpha in the line search (> 0) 
%
% - sfgrd (real scalar, optional, default value 0.2): safeguard parameter
%   for the line search. to avoid numerical problems that can occur with
%   the quadratic interpolation if the derivative at one endpoint is too
%   large w.r.t. the one at the other (which leads to choosing a point
%   extremely near to the other endpoint), a *safeguarded* version of
%   interpolation is used whereby the new point is chosen in the interval
%   [ as * ( 1 + sfgrd ) , am * ( 1 - sfgrd ) ], being [ as , am ] the
%   current interval, whatever quadratic interpolation says. If you
%   experiemce problems with the line search taking too many iterations to
%   converge at "nasty" points, try to increase this
%
% - mina (real scalar, optional, default value 1e-16): if the algorithm
%   determines a stepsize value <= mina, this is taken as an indication
%   that something has gone wrong (the gradient is not a direction of
%   descent, so maybe the function is not differentiable) and computation
%   is stopped. It is legal to take mina = 0, thereby in fact skipping this
%   test.
%
% Output:
%
% - p (real scalar): the best (largest) value of the dual function found
%   so far (possibly the optimal one), i.e., the best available lower
%   bound on the optimal value of (P)
%
% - status (string, optional): a string describing the status of the
%   algorithm at termination, with the following possible values:
%
%   = 'optimal': the algorithm terminated having proven that x is a(n
%     approximately) optimal solution, i.e., the norm of the gradient at x
%     is less than the required threshold
%
%   = 'stopped': the algorithm terminated having exhausted the maximum
%     number of iterations: x is the bast solution found so far, but not
%     necessarily the optimal one
%
%   = 'error': the algorithm found a numerical error that prevents it from
%     continuing optimization (see mina above)
%
% - l ([ 2 * n x 1 ] real column vector, optional): the best Lagrangian
%   multipliers found so far (possibly the optimal ones)
%
% - v (real scalar): the best (smallest) value of the objective function
%   found so far, that of the best feasible solution found so far by the
%   Lagrangian heuristic (see x below), i.e., the best available upper
%   bound on the optimal value of (P); if dolh == false, then v == Inf
%
% - x ([ n x 1 ] real column vector, optional): the best feasible solution
%   found so far by the Lagrangian heuristic (possibly the optimal one);
%   if dolh == false, then x == []
%
%{
 =======================================
 Author: Antonio Frangioni
 Date: 09-09-21
 Version 1.10
 Copyright Antonio Frangioni
 =======================================
%}

Interactive = false;  % if we pause at every iteration

Verbose = true;      % if we log every iteration

% reading and checking input- - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if ~ isstruct( BCQP )
   error( 'BCQP not a struct' );
end

if ~ isfield( BCQP , 'Q' ) || ~ isfield( BCQP , 'q' ) || ...
   ~ isfield( BCQP , 'u' )
   error( 'BCQP not a well-formed struct' );
end

if ~ isreal( BCQP.Q ) || ~ isreal( BCQP.q ) || ~ isreal( BCQP.u )
   error( 'BCQP not a well-formed struct' );
end

n = size( BCQP.q , 1 );
if size( BCQP.q , 2 ) ~= 1 || ~ isequal( size( BCQP.Q ) , [ n n ] ) || ...
   ~ isequal( size( BCQP.u ) , [ n 1 ] )
   error( 'BCQP not a well-formed struct' );
end

% compute straight away the Cholesky factorization of Q, this will be used
% at each iteration to solve the Lagrangian relaxation
if ~ Interactive  % and since this is costly, include it in the time
   tic;
end

[ R , h ] = chol( BCQP.Q );
if h > 0
   error( 'BCQP.Q not positive definite, this is not supported (yet)' );
end

if ~isempty( varargin )
   dolh = logical( varargin{ 1 } );
else
   dolh = true;
end

if length( varargin ) > 1
   eps = varargin{ 2 };
   if ~ isreal( eps ) || ~ isscalar( eps )
      error( 'eps is not a real scalar' );
   end
else
   eps = 1e-6;
end

if length( varargin ) > 2
   MaxFeval = round( varargin{ 3 } );
   if ~ isscalar( MaxFeval )
      error( 'MaxFeval is not an integer scalar' );
   end
else
   MaxFeval = 1000;
end

if length( varargin ) > 3
   m1 = varargin{ 4 };
   if ~ isscalar( m1 )
      error( 'm1 is not a real scalar' );
   end
   if m1 <= 0 || m1 >= 1
      error( 'm1 is not in (0 ,1)' );
   end       
else
   m1 = 0.01;
end

if length( varargin ) > 4
   m2 = varargin{ 5 };
   if ~ isscalar( m1 )
      error( 'm2 is not a real scalar' );
   end
   if m2 <= 0 || m2 >= 1
      error( 'm2 is not in (0, 1)' );
   end
else
   m2 = 0.9;
end

if length( varargin ) > 5
   astart = varargin{ 6 };
   if ~ isscalar( astart )
      error( 'astart is not a real scalar' );
   end
   if astart < 0
      error( 'astart must be > 0' );
   end       
else
   astart = 1;
end

if length( varargin ) > 6
   sfgrd = varargin{ 7 };
   if ~ isscalar( sfgrd )
      error( 'sfgrd is not a real scalar' );
   end
   if sfgrd <= 0 || sfgrd >= 1
      error( 'sfgrd is not in (0, 1)' );
   end
else
   sfgrd = 0.2;
end

if length( varargin ) > 7
   mina = varargin{ 8 };
   if ~ isscalar( mina )
      error( 'mina is not a real scalar' );
   end
   if mina < 0
      error( 'mina is < 0' );
   end       
else
   mina = 1e-12;
end

% initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

fprintf( 'Lagrangian Dual\n');
if Verbose
   if dolh
      fprintf( 'feval\tub\t\tp(l)\t\tgap\t\tls feval\ta*\n\n' );
   else
      fprintf( 'feval\tp(l)\t\t|| grad ||\tls feval\ta*\n\n' );
   end
end

x = [];
v = Inf;
feval = 0;
i = 1;

lambda = zeros( 2 * n , 1 );
[ p , lastg ] = phi( lambda );

% main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

while true
    % project the direction = - gradient over the active constraints- - - -
    d = - lastg;
    d( lambda <= 1e-12 & d < 0 ) = 0;

    % output statistics - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
    if dolh
       % compute the relative gap
       gap = ( v + p ) / max( [ abs( v ) 1 ] );

       if Verbose
          fprintf( '%4d\t%1.8e\t%1.8e\t%1.4e\t' , feval , v , - p , gap );
       end

       if gap <= eps   % look, ma! the *true* stopping criterion, for once!
          status = 'optimal';
          break;
       end
    else
       % compute the norm of the projected gradient
       gnorm = norm( d );

       if Verbose
          fprintf( '%4d\t%1.8e\t%1.4e\t' , feval , - p , gnorm );
       end

       if feval == 1
          gnorm0 = gnorm;
       end
       if gnorm <= eps * gnorm0
          status = 'optimal';
          break;
       end
    end
    % stopping criteria - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    if feval > MaxFeval
       status = 'stopped';
       break;
    end

    % compute step size - - - - - - - - - - - - - - - - - - - - - - - - - -

    % first, compute the maximum feasible stepsize maxt such that
    %
    %   0 <= lambda( i ) + maxt * d( i )   for all i

    ind = d < 0;  % negative gradient entries
    if any( ind )
       maxt = min( astart , min( - lambda( ind ) ./ d( ind ) ) );
    else
       maxt = astart;
    end

    % now run the line search
    phip0 = lastg' * d;    
    [ a , p ] = ArmijoWolfeLS( p , phip0 , maxt , m1 , m2 );

    if Verbose
       fprintf( '\t%1.4e\n' , a );
    end
    if a <= mina
       fprintf( '\tERR\n' );
       status = 'error';
       break;
    end
    
    % compute new point - - - - - - - - - - - - - - - - - - - - - - - - - -

    lambda = lambda + a * d;

    % iterate - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    i = i + 1;

    if Interactive
       pause;
    end
end

% end of main loop- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

p = -p;  % the value of the dual function is changed sign inside

if dolh
   fprintf( 'stop: %d iter, status = %s, lbest = %1.8e, gap = %1.4e\n' ,...
            i , status , p , gap );
else
   fprintf( 'stop: %d iter, status = %s, lbest = %1.8e\n' , ...
            i , status , p );
end

if ~ Interactive
   toc
end

if nargout > 1
   varargout{ 1 } = status;
end

if nargout > 2
   varargout{ 2 } = lambda;
end

if nargout > 3
   varargout{ 3 } = v;
end

if nargout > 4
   varargout{ 4 } = x;
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% inner functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function [ p , y ] = solveLagrangian( lambda )
% The Lagrangian relaxation of the problem is
%
%  min { (1/2) x' * Q * x + q * x - lambda^+ * ( u - x ) - lambda^- * ( x )
% = min { (1/2) x' * Q * x + ( q + lambda^+ - lambda^- ) * x - lambda^+ * u
%
% where lambda^+ are the first n components of lambda[], and lambda^- the
% last n components; both are constrained to be >= 0
%
% The optimal solution of the Lagrangian relaxation is the (unique)
% solution of the linear system
%
%       Q * x = - q - lambda^+ + lambda^-
%
% Since we have computed at the beginning the Cholesky factorization of Q,
% i.e., Q = R' * R, where R is upper triangular and therefore RT is lower
% triangular, we obtain this by just two triangular backsolves:
%
%       R' * z = - q - lambda^+ + lambda^-
%
%       R * x = z
%
% return the function value and the primal solution

ql = BCQP.q + lambda( 1 : n ) - lambda( n + 1 : end );
opts.LT = true;
z = linsolve( R' , - ql , opts );
opts.LT = false;
opts.UT = true;
y = linsolve( R , z , opts );

% compute phi-value
p = ( 0.5 * y' * BCQP.Q + ql' ) * y - lambda( 1 : n )' * BCQP.u;

feval = feval + 1;

end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function [ p , g ] = phi( lambda )
% phi( lambda ) is the Lagrangian function of the problem. With x the
% optimal solution of the minimization problem (see solveLagrangian()), the
% gradient at lambda is [ x - u ; - x ]
%
% however, the line search is written for minimization but we rather want
% to maximize phi(), hence we have to change the sign of both function
% values and gradient entries

% solve the Lagrangian relaxation
[ p , y ] = solveLagrangian( lambda );
p = - p;

% compute gradient
g = [ BCQP.u - y ; y ];

if dolh
   % compute an heuristic solution out of the solution y of the Lagrangian
   % relaxation by projecting y on the box

   y( y < 0 ) = 0;
   ind = y > BCQP.u;
   y( ind ) = BCQP.u( ind );

   % compute cost of feasible solution
   pv = 0.5 * y' * BCQP.Q * y + BCQP.q' * y;

   if pv < v   % it is better than best one found so far
      x = y;   % y becomes the incumbent
      v = pv;
   end
end
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function [ p , pp ] = phi1d( alpha )
% phi( lambda ) is the Lagrangian function of the problem; then,
% phi( alpha ) = phi( lambda + alpha d ) and
% phi'( alpha ) = < \nabla phi( lambda + alpha * d ) , d >

[ p , lastg ] = phi( lambda + alpha * d );
pp = d' * lastg;

end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function [ a , phia ] = ArmijoWolfeLS( phi0 , phip0 , as , m1 , m2 )

% performs an Armijo-Wolfe Line Search.
%
% phi0 = phi( 0 ), phip0 = phi'( 0 ) < 0
%
% as > 0 is the first value to be tested, and it is also the *maximum*
% possible stepsize: if phi'( as ) < 0 then the LS is immediately
% terminated
%
% m1 and m2 are the standard Armijo-Wolfe parameters; note that the strong
% Wolfe condition is used
%
% returns the optimal step and the optimal f-value

a = as;
[ phia , phips ] = phi1d( a );
if phips <= 0
   if Verbose
      fprintf( '%2d' , 1 );
   end
   return;
end

lsiter = 1;  % count ls iterations

am = 0;
phipm = phip0;
while ( feval <= MaxFeval ) && ( ( as - am ) ) > mina && ( phips > 1e-12 )

   % compute the new value by safeguarded quadratic interpolation
   a = ( am * phips - as * phipm ) / ( phips - phipm );
   a = max( [ am * ( 1 + sfgrd ) min( [ as * ( 1 - sfgrd ) a ] ) ] );

   % compute phi( a )
   [ phia , phip ] = phi1d( a );

   if ( phia <= phi0 + m1 * a * phip0 ) && ( abs( phip ) <= - m2 * phip0 )
      break;  % Armijo + strong Wolfe satisfied, we are done
   end

   % restrict the interval based on sign of the derivative in a
   if phip < 0
      am = a;
      phipm = phip;
   else
      as = a;
      if as <= mina
         break;
      end
      phips = phip;
   end
   lsiter = lsiter + 1;
end

if Verbose
   fprintf( '%2d' , lsiter );
end

end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end  % the end- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




