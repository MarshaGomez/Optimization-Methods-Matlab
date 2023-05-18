function [ v , varargout ] =  PGBCQP( BCQP , varargin )

%function [ v , x , status ] = PGBCQP( BCQP , walg , alpha , eps ,
%                                      MaxIter )
%
% Apply the Projected Gradient algorithm with exact line search to the
% convex Box-Constrained Quadratic program
%
%  (P) min { (1/2) x^T * Q * x + q * x : 0 <= x <= u }
%
% encoded in the structure BCQP.
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
% - walg (single character, optional, default value 's') if walg == 's' or
%   'S' then the standard version of the projected gradient method is used
%   that first projects the gradient on the feasible direction cone and
%   then choses the stepsize along the projected gradient. if, rather,
%   walg == 'g' or 'G', then Goldstein's version of the projected gradient
%   method is used that first choses the stepsize along the un-projected
%   gradient and then projects the iterate back onto the feasible region.
%
% - alpha (real scalar, optional, default value 0): if > 0 then it is the
%   Fixed Stepsize to be used along the chosen descent direction (whatever
%   of the two it is), otherwise an exact LS (this being a quadratic
%   objective) is used to select the stepsize.
%
% - eps (real scalar, optional, default value 1e-6): the accuracy in the 
%   stopping criterion: the algorithm is stopped when the norm of the
%   (projected) gradient is less than or equal to eps
%
% - MaxIter (integer scalar, optional, default value 1000): the maximum
%   number of iterations
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
%     number of iterations: x is the bast solution found so far, but not
%     necessarily the optimal one
%
%
%{
 =======================================
 Author: Antonio Frangioni
 Date: 08-09-21
 Version 2.00
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

if ~isempty( varargin )
   walg = varargin{ 1 };
   if ~ ischar( walg )
      error( 'walg is not a character' );
   end
   walg = lower( walg );
   if walg == 's'
      Goldstein = false;
   else
      Goldstein = true;
   end
else
   Goldstein = false;
end

if length( varargin ) > 1
   alpha = varargin{ 2 };
   if ~ isreal( alpha ) || ~ isscalar( alpha )
      error( 'alpha is not a real scalar' );
   end
else
   alpha = 0;
end

if length( varargin ) > 2
   eps = varargin{ 3 };
   if ~ isreal( eps ) || ~ isscalar( eps )
      error( 'eps is not a real scalar' );
   end
else
   eps = 1e-6;
end

if length( varargin ) > 3
   MaxIter = round( varargin{ 4 } );
   if ~ isscalar( MaxIter )
      error( 'MaxIter is not an integer scalar' );
   end
else
   MaxIter = 1000;
end

% initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

x = BCQP.u / 2;  % start from the middle of the box

fprintf( 'Projected gradient method, ' );
if Goldstein
   fprintf( 'Goldstein''s version\n' );
else
   fprintf( 'standard version\n' );
end
if Verbose
   fprintf( 'iter\tf(x)\t\t\t||nabla f(x)||\n\n' );
end

i = 1;

% main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if ~ Interactive
   tic;
end

while true

    % compute gradient- - - - - - - - - - - - - - - - - - - - - - - - - - -

    g = BCQP.Q * x + BCQP.q;

    if Goldstein
       if alpha <= 0
          % compute optimal unbounded stepsize:
          % a = ( g' * g ) / ( g' * Q * g )

          den = g' * BCQP.Q * g;
    
          if den <= 1e-16  % d' * Q * d = 0  ==>  f is linear along d
             t = 1;     % just take any stepsize (1)
          else
             % optimal unbounded stepsize restricted to max feasible step
             t = ( g' * g ) / den;
          end
       else
          t = alpha;  % fixed stepsize
       end
  
       y = x - t * g;  % do the step to a possibly unfeasible point
 
       % project the possibly unfeasible point back into the box
       y = max( 0 , min( BCQP.u , y ) );

       % compute the norm of the movement
       ng = norm( y - x );

       % actually make the step
       x = y;
    else
       d = - g;
  
       % project the direction over the active constraints
       d( BCQP.u - x <= 1e-12 & d > 0 ) = 0;
       d( x <= 1e-12 & d < 0 ) = 0;

       % compute the norm of the (projected) gradient
       ng = norm( d );
    end

    % output statistics - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if Verbose
       v = 0.5 * x' * BCQP.Q * x + BCQP.q' * x;
       fprintf( '%4d\t%1.8e\t\t%1.4e\n' , i , v , ng );
    end
   
    % stopping criteria - - - - - - - - - - - - - - - - - - - - - - - - - -
    if ng <= eps
       status = 'optimal';
       break;
    end
    
    if i > MaxIter
       status = 'stopped';
       break;
    end

    if ~ Goldstein
       % compute step size- - - - - - - - - - - - - - - - - - - - - - - - -

       if alpha <= 0
          % first, compute the maximum feasible stepsize maxt such that
          %
          %   0 <= x( i ) + maxt * d( i ) <= u( i )   for all i

          ind = d > 0;  % positive gradient entries
          maxt = min( ( BCQP.u( ind ) - x( ind ) ) ./ d( ind ) );
          ind = d < 0;  % negative gradient entries
          maxt = min( [ maxt min( - x( ind ) ./ d( ind ) ) ] );

          % compute optimal unbounded stepsize:
          % min (1/2) ( x + t d )' * Q * ( x + t d ) + q' * ( x + t d ) =
          %     (1/2) t^2 ( d' * Q * d ) + t d' * ( Q * x + q ) [ + c ]
          %
          % ==> t = - d' * ( Q * x + q ) / d' * Q * d
          %       = - g' * d / d' * Q * d 
          %     (note that - g' * d < 0, hence t > 0 as expected)
          %
          den = d' * BCQP.Q * d;
    
          if den <= 1e-16  % d' * Q * d = 0  ==>  f is linear along d
             t = maxt;     % just take the maximum possible stepsize
          else
             % optimal unbounded stepsize restricted to max feasible step
             t = min( [ ( - g' * d ) / den , maxt ] );
          end
       else
          t = alpha;  % fixed stepsize
       end

       % compute new point- - - - - - - - - - - - - - - - - - - - - - - - -

       x = x + t * d;
    end

    % iterate - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    i = i + 1;

    if Interactive
       pause;
    end
end

% end of main loop- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if ~ Verbose
   v = 0.5 * x' * BCQP.Q * x + BCQP.q' * x;
end

fprintf( 'stop: %d iter, status = %s, fbest = %1.8e\n' , ...
         i , status , v );
if ~ Interactive
   toc
end

if nargogolut > 1
   varargout{ 1 } = x;
end

if nargout > 2
   varargout{ 2 } = status;
end

end  % the end- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




