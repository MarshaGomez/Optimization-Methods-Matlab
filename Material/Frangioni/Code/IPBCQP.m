function [ v , varargout ] =  IPBCQP( BCQP , varargin )

%function [ v , x , status ] = IPBCQP( BCQP , iD , MaxIter , eps )
%
% Apply the Primal-Dual (feasible) Interior (barrier) Method to the convex
% Box-Constrained Quadratic program
%
%  (P) min { (1/2) x^T * Q * x + q * x : 0 <= x <= u }
%
% encoded in the structure BCQP, where Q must be strictly positive
% definite.
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
% - id (real scalar, optional, default value 1): the initial dual solution
%   is chosen as a scalar multiple of u by some value t. t is chosen as
%   the maximum value of u (so that it is "well scaled" w.r.t. the data of
%   the problem) multiplied by id. A large id means we start "a lot inside"
%   the dual polyhedron which may give a bad dual solution, but a too small
%   id means that we start "too close to the boundary" of the dual
%   polyhedron and this may impair convergence speed.
%
% - MaxIter (integer scalar, optional, default value 100): the maximum
%   number of iterations
%
% - eps (real scalar, optional, default value 1e-10): the accuracy in the 
%   stopping criterion: the algorithm is stopped when the relative gap
%   between the value of the current primal and dual feasible solutions is
%   less than or equal to eps
%
% Output:
%
% - v (real scalar): the best function value found so far (possibly the
%   optimal one)
%
% - x ([ n x 1 ] real column vector, optional): the best solution found so
%   far (possibly the optimal one)
%
% - status (string, optional): a string describing the status of the
%   algorithm at termination, with the following possible values:
%
%   = 'optimal': the algorithm terminated having proven that x is a(n
%     approximately) optimal solution, i.e., the norm of the gradient at x
%     is less than the required threshold
%
%   = 'stopped': the algorithm terminated having exhausted the maximum
%     number of iterations: x is the best solution found so far, but not
%     necessarily the optimal one
%{
 =======================================
 Author: Antonio Frangioni
 Date: 03-12-21
 Version 1.20
 Copyright Antonio Frangioni
 =======================================
%}

printCS = 0;

Interactive = false; % if we pause at every iteration

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

if ~isempty( varargin )
   id = varargin{ 1 };
   if ~ isscalar( id )
      error( 'id is not a real scalar' );
   end
   if id <= 0
      error( 'id is not > 0' );
   end
else
   id = 1;
end

if length( varargin ) > 1
   MaxIter = round( varargin{ 2 } );
   if ~ isscalar( MaxIter )
      error( 'MaxIter is not an integer scalar' );
   end
else
   MaxIter = 1000;
end

if length( varargin ) > 2
   eps = varargin{ 3 };
   if ~ isreal( eps ) || ~ isscalar( eps )
      error( 'eps is not a real scalar' );
   end
else
   eps = 1e-10;
end

% initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% the Slackened KKT System for (P) (written without slacks) is
%
%   Q x + q + \lambda^+ - \lambda^- = 0
%
%   \lambda^+ ( u - x ) = \mu * e
%
%   \lambda^-     x     = \mu * e
%
%   0 <= x <= u
%
%   \lambda^+ >= 0   ,   \lambda^- >= 0
%
% where e is the vector of all ones.
%
% if x and ( \lambda^+ , \lambda^- ) satisfy SKKTS, then
%
%   v = (1/2) x' Q x + q x
%
%   p = - \lambda^+ u - (1/2) x' Q x
%
% are, respectively, a valid upper and lower bound on v(P), and
%
%   v - p = \mu * 2 * n
%
% With x -> x + dx, \lambda^+ -> lp + dlp, \lambda^+ -> lm + dlm such that
%
%   Q x + q + lp - lm = 0
%
% SKKTS becomes
%
%   Q dx + dlp - dlm = 0                             (1)
%
%   ( lp + dlp ) ( u - x - dx ) = \mu * e            (2)
%
%   ( lm + dlm ) (   x   + dx ) = \mu * e            (3)
%
%   - x <= dx <= u - x                               (4)
%
%   dlp >= - lp   ,   dlm >= - lp                    (5)
%
% inequalities (4) and (5) are just the bounds, and will be taken care of
% by the appropriate choice of the stepsize. the well-known trick is
% linearizing the nonlinear equalities (2) and (3) by just ignoring the
% bilinear terms, which leads to
%
%   Q dx + dlp - dlm = 0                             (1)
%
% - lp dx + dlp ( u - x ) = \mu * e - lp ( u - x )   (2')
%
%   lm dx + dlm     x     = \mu * e - lm     x       (3')
%
% we can then use (2') and (3') to write
%
%    dlp = ( \mu * e + lp dx ) ./ ( u - x ) - lp     (2'')
%
%    dlm = ( \mu * e - lm dx ) ./     x     - lm     (3'')
%
% putting (2'') and (3'') in (1) gives
%
%   H dx = w                                         (1')
%
% with 
%
%   H = Q + lp ./ ( u - x ) + lm ./ x
%
%   w = \mu * [ e ./ ( u - x ) - e ./ x ] + lp - lm
%
% where note that lp - lm = - Q x - q
%
% The term H - Q is diagonal and strictly positive, hence H is strictly
% positive definite and nonsingular and (1') has a unique solution.
%
% To initialize the algorithm we take x straight in the middle of the box,
% and then it would be simple to satisfy
%
%   Q x + q + lp - lm = 0
%
% by taking lm = [ Q x + q ]_+ and lp = [ - Q x - q ]_+. however, by doing
% so lm and lp would not be interior. The obvious solution is to add to
% both a term eps * e with some small eps (1e-6)

% compute a feasible interior primal solution (the middle of the box)
x = BCQP.u / 2;

% compute a feasible interior dual solution satifying SKKTS with x for
% some \mu we don't care much of
g = BCQP.Q * x + BCQP.q;

% initialise lp and lm as a scalar multiple of e
lp = id * max( BCQP.u ) * ones( n , 1 );  
lm = lp;
ind = g >= 0;
lm( ind ) = lm( ind ) + g( ind );
ind = ~ind;
lp( ind ) = lp( ind ) - g( ind );

fprintf( 'Primal-Dual Interior Point method\n');
if Verbose
   fprintf( 'iter\tv\t\tp\t\tgap\n\n' );
end

i = 1;

% main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if ~ Interactive
   tic;
end

while true

    % compute upper and lower bound - - - - - - - - - - - - - - - - - - - -

    xQx = x' * BCQP.Q * x;
    v = 0.5 * xQx + BCQP.q' * x;
    p = - lp' * BCQP.u - 0.5 * xQx;
    gap = ( v - p ) / max( abs( v ) , 1 );

    % output statistics - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if Verbose
       fprintf( '%4d\t%1.8e\t%1.8e\t%1.4e\n' , i , v , p , gap );
   
       if printCS
          fprintf( 'x\tlm\tu-x\tlp\n' );
          for h = 1 : n
              fprintf( '%1.1e\t%1.1e\t%1.1e\t%1.1e\n' , x( h ) , ...
                       lm( h ) , BCQP.u( h ) - x( h ) , lp( h ) );
          end
       end
    end

    % stopping criteria - - - - - - - - - - - - - - - - - - - - - - - - - -

    if gap <= eps      % look, ma! the *true* stopping criterion, for once!
       status = 'optimal';
       break;
    end
    
    if i > MaxIter
       status = 'stopped';
       break;
    end

    % solve the SKKTS - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % note: the "complicated" term in W has the form
    %
    %  mu [ 1 / ( u_i - x_i ) - 1 / x_i ]
    %
    % which can be rewritten
    %
    %  mu ( u_i - 2 x_i ) / [ ( u_i - x_i ) * x_i ]
    %
    % it appears this last form is *vastly* more numerically stable, try
    % for yourself if you don't believe me

    mu = ( v - p ) / ( 4 * n * n );  % use \rho = 1 / ( # of constraints )

    umx = BCQP.u - x;
    H = BCQP.Q + diag( lp ./ umx + lm ./ x );
    % w = mu * ( ones( n , 1 ) ./ umx - ones( n , 1 ) ./ x ) + lp - lm;
    w = mu * ( BCQP.u - 2 * x ) ./ ( umx .* x ) + lp - lm;

    opts.SYM = true;     % tell it H is positive definite (hence of course
    opts.POSDEF = true;  % symmetric), it'll probably do the right thing
                         % and use Cholesky to solve the system
    dx = linsolve( H , w , opts );

    dlp = ( mu * ones( n , 1 ) + lp .* dx ) ./ umx - lp;

    dlm = ( mu * ones( n , 1 ) - lm .* dx ) ./  x  - lm;    
    
    % compute maximum feasible primal stepsize- - - - - - - - - - - - - - -

    ind = dx < 0;  % negative direction entries
    if any( ind )
       maxt = min( - x( ind ) ./ dx( ind ) );
    else
       maxt = inf;
    end
    ind = dx > 0;  % positive direction entries
    if any( ind )
       maxt = min( maxt , min( umx( ind ) ./ dx( ind ) ) );
    end
    
    % compute maximum feasible dual stepsize- - - - - - - - - - - - - - - -

    ind = dlp < 0;  % negative direction entries
    if any( ind )
       maxt = min( maxt , min( - lp( ind ) ./ dlp( ind ) ) );
    end
    ind = dlm < 0;  % negative direction entries
    if any( ind )
       maxt = min( maxt , min( - lm( ind ) ./ dlm( ind ) ) );
    end

    % compute new primal-dual solution- - - - - - - - - - - - - - - - - - -

    maxt = 0.9995 * maxt;  % ensure the new iterate remains interior

    x = x + maxt * dx;
    lp = lp + maxt * dlp;
    lm = lm + maxt * dlm;

    % iterate - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    i = i + 1;

    if Interactive
       pause;
    end
end

% end of main loop- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

fprintf( 'stop: %d iter, status = %s, fbest = %1.8e\n' , ...
         i , status , v );
if ~ Interactive
   toc
end

if nargout > 1
   varargout{ 1 } = x;
end

if nargout > 2
   varargout{ 2 } = status;
end

end  % the end- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




