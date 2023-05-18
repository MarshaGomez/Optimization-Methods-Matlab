function [ x , status ] =  NCG( f , varargin )

%function [ x , status ] = NCG( f , x , wf , rstart , eps , astart , 
%                               MaxFeval , m1 , m2 , tau , sfgrd , MInf ,
%                               mina )
%
% Apply a Nonlinear Conjugated Gradient algorithm for the minimiztion of
% the provided function f, which must have the following interface:
%
%   [ v , g ] = f( x )
%
% Input:
%
% - x is either a [ n x 1 ] real (column) vector denoting the input of
%   f(), or [] (empty).
%
% Output:
%
% - v (real, scalar): if x == [] this is the best known lower bound on
%   the unconstrained global optimum of f(); it can be -Inf if either f()
%   is not bounded below, or no such information is available. If x ~= []
%   then v = f(x).
%
% - g (real, [ n x 1 ] real vector): this also depends on x. if x == []
%   this is the standard starting point from which the algorithm should
%   start, otherwise it is the gradient of f() at x (or a subgradient if
%   f() is not differentiable at x, which it should not be if you are
%   applying the gradient method to it).
%
% The other [optional] input parameters are:
%
% - x (either [ n x 1 ] real vector or [], default []): starting point.
%   If x == [], the default starting point provided by f() is used.
%
% - wf (integer scalar, optional, default value 0): which of the Nonlinear
%   Conjugated Gradient formulae to use. Possible values are:
%   = 0: Fletcher-Reeves
%   = 1: Polak-Ribiere
%   = 2: Hestenes-Stiefel
%   = 3: Dai-Yuan
%
% - rstart (integer scalar, optional, default value 0): if > 0, restarts
%   (setting beta = 0) are performed every n * rstart iterations
%
% - eps (real scalar, optional, default value 1e-6): the accuracy in the 
%   stopping criterion: the algorithm is stopped when the norm of the
%   gradient is less than or equal to eps. If a negative value is provided,
%   this is used in a *relative* stopping criterion: the algorithm is
%   stopped when the norm of the gradient is less than or equal to
%   (- eps) * || norm of the first gradient ||.
%
% - astart (real scalar, optional, default value 1): starting value of
%   alpha in the line search (> 0) 
%
% - MaxFeval (integer scalar, optional, default value 1000): the maximum
%   number of function evaluations (hence, iterations will be not more than
%   MaxFeval because at each iteration at least a function evaluation is
%   performed, possibly more due to the line search).
%
% - m1 (real scalar, optional, default value 0.01): first parameter of the
%   Armijo-Wolfe-type line search (sufficient decrease). Has to be in (0,1)
%
% - m2 (real scalar, optional, default value 0.9): typically the second
%   parameter of the Armijo-Wolfe-type line search (strong curvature
%   condition). It should to be in (0,1); if not, it is taken to mean that
%   the simpler Backtracking line search should be used instead
%
% - tau (real scalar, optional, default value 0.0): scaling parameter for
%   the line search. In the Armijo-Wolfe line search it is used in the
%   first phase: if the derivative is not positive, then the step is
%   divided by tau (which is < 1, hence it is increased). In the 
%   Backtracking line search, each time the step is multiplied by tau
%   (hence it is decreased).
%
% - sfgrd (real scalar, optional, default value 0.01): safeguard parameter
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
% - MInf (real scalar, optional, default value -Inf): if the algorithm
%   determines a value for f() <= MInf this is taken as an indication that
%   the problem is unbounded below and computation is stopped
%   (a "finite -Inf").
%
% - mina (real scalar, optional, default value -1e-16): if the algorithm
%   determines a stepsize value <= abs( mina ), this is taken as an
%   indication that something has gone wrong (the direction is not of
%   descent), and some drastic action need be taken. Which action is
%   decided by the sign of mina: if mina >= 0 then the computation is
%   stopped, otherwise (by default) a reset of the algorithm is performed,
%   equivalent to forcing beta = 0 at the next iteration
%
% Output:
%
% - x ([ n x 1 ] real column vector): the best solution found so far.
%
% - status (string): a string describing the status of the algorithm at
%   termination 
%
%   = 'optimal': the algorithm terminated having proven that x is a(n
%     approximately) optimal solution, i.e., the norm of the gradient at x
%     is less than the required threshold
%
%   = 'unbounded': the algorithm has determined an extrenely large negative
%     value for f() that is taken as an indication that the problem is
%     unbounded below (a "finite -Inf", see MInf above)
%
%   = 'stopped': the algorithm terminated having exhausted the maximum
%     number of iterations: x is the bast solution found so far, but not
%     necessarily the optimal one
%
%   = 'error': the algorithm found a numerical error that prevents it from
%     continuing optimization (see mina above)
%
%{
 =======================================
 Author: Antonio Frangioni
 Date: 12-11-22
 Version 1.40
 Copyright Antonio Frangioni
 =======================================
%}

Plotf = 1;
% 0 = nothing is plotted
% 1 = the level sets of f and the trajectory are plotted (when n = 2)
% 2 = the function value / gap are plotted, iteration-wise
% 3 = the function value / gap are plotted, function-evaluation-wise

Interactive = false;  % if we pause at every iteration

% reading and checking input- - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if ~ isa( f , 'function_handle' )
   error( 'f not a function' );
end

if isempty( varargin ) || isempty( varargin{ 1 } )
   [ fStar , x ] = f( [] );
else
   x = varargin{ 1 };
   if ~ isreal( x )
      error( 'x not a real vector' );
   end

   if size( x , 2 ) ~= 1
      error( 'x is not a (column) vector' );
   end
   
   fStar = f( [] );
end

n = size( x , 1 );

if length( varargin ) > 1
   wf = round( varargin{ 2 } );
   if ~ isscalar( wf )
      error( 'wf is not a real scalar' );
   end
   if wf < 0 || wf > 3
      error( 'unknown NCG formula %d' , wf );
   end   
else
   wf = 0;
end

if length( varargin ) > 2
   rstart = round( varargin{ 3 } );
   if ~ isscalar( rstart )
      error( 'rstart is not an integer scalar' );
   end
else
   rstart = 0;
end

if length( varargin ) > 3
   eps = varargin{ 4 };
   if ~ isreal( eps ) || ~ isscalar( eps )
      error( 'eps is not a real scalar' );
   end
else
   eps = 1e-6;
end

if length( varargin ) > 4
   astart = varargin{ 5 };
   if ~ isscalar( astart )
      error( 'astart is not a real scalar' );
   end
   if astart < 0
      error( 'astart must be > 0' );
   end       
else
   astart = 1;
end

if length( varargin ) > 5
   MaxFeval = round( varargin{ 6 } );
   if ~ isscalar( MaxFeval )
      error( 'MaxFeval is not an integer scalar' );
   end
else
   MaxFeval = 1000;
end

if length( varargin ) > 6
   m1 = varargin{ 7 };
   if ~ isscalar( m1 )
      error( 'm1 is not a real scalar' );
   end
   if m1 <= 0 || m1 >= 1
      error( 'm1 is not in (0 ,1)' );
   end       
else
   m1 = 0.01;
end

if length( varargin ) > 7
   m2 = varargin{ 8 };
   if ~ isscalar( m1 )
      error( 'm2 is not a real scalar' );
   end
else
   m2 = 0.9;
end

AWLS = ( m2 > 0 && m2 < 1 );

if length( varargin ) > 8
   tau = varargin{ 9 };
   if ~ isscalar( tau )
      error( 'tau is not a real scalar' );
   end
   if tau <= 0 || tau >= 1
      error( 'tau is not in (0 ,1)' );
   end       
else
   tau = 0.9;
end

if length( varargin ) > 9
   sfgrd = varargin{ 10 };
   if ~ isscalar( sfgrd )
      error( 'sfgrd is not a real scalar' );
   end
   if sfgrd <= 0 || sfgrd >= 1
      error( 'sfgrd is not in (0, 1)' );
   end
else
   sfgrd = 0.01;
end

if length( varargin ) > 10
   MInf = varargin{ 11 };
   if ~ isscalar( MInf )
      error( 'MInf is not a real scalar' );
   end
else
   MInf = - Inf;
end

if length( varargin ) > 11
   mina = varargin{ 12 };
   if ~ isscalar( mina )
      error( 'mina is not a real scalar' );
   end
else
   mina = - 1e-16;
end

% "global" variables- - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

lastx = zeros( n , 1 );  % last point visited in the line search
lastg = zeros( n , 1 );  % gradient of lastx
d = zeros( n , 1 );      % NGC's direction
feval = 0;               % f() evaluations count ("common" with LSs)

% initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

fprintf( 'NCG method ' );
switch wf
   case 0
      fprintf( '(Fletcher-Reeves)\n' );
   case 1
      fprintf( '(Polak-Ribiere)\n' );
   case 2
      fprintf( '(Hestenes-Stiefel)\n' );
   otherwise
      fprintf( '(Dai-Yuan)\n' );
end
if fStar > - Inf
   fprintf( 'feval\trel gap');
else
   fprintf( 'feval\tf(x)');
end
fprintf( '\t    || g(x) ||\tbeta\t ls feval   a*\n\n');

if Plotf > 1
   % expected number of iterations for the gap plot
   if Plotf == 2
      MaxIter = MaxFeval / 3;
      ylim( [ 1e-15 1e+1 ] );
   else
      MaxIter = MaxFeval;
      ylim( [ 1e-15 1e+4 ] );
    end
   gap = [];

   xlim( [ 0 MaxIter ] );
   ax = gca;
   ax.FontSize = 16;
   ax.Position = [ 0.03 0.07 0.95 0.92 ];
   ax.Toolbar.Visible = 'off';
end

v = f2phi( 0 );
ng = norm( lastg );
if eps < 0
   ng0 = - ng;  % norm of first subgradient: why is there a "-"? ;-)
else
   ng0 = 1;     % un-scaled stopping criterion
end

sgiter = true;  % true if there is no previous direction, i.e., a gradient
                % step will be taken (beta will not even be computed)
iter = 1;       % iterations, as distinct from function evaluations

% main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

while true
    % output statistics - - - - - - - - - - - - - - - - - - - - - - - - - -
  
    if fStar > - Inf
       gapk = ( v - fStar ) / max( [ abs( fStar ) 1 ] );

       fprintf( '%4d\t%1.4e  %1.4e' , feval , gapk , ng );

       if Plotf > 1
          if Plotf == 2
             gap( end + 1 ) = gapk;
          end
          semilogy( gap , 'Color' , 'k' , 'LineWidth' , 2 );
          drawnow;
       end
    else
       fprintf( '%4d\t%1.4e\t  %1.4e' , feval , v , ng );

       if Plotf > 1
          if Plotf == 2
             gap( end + 1 ) = v;
          end
          plot( gap , 'Color' , 'k' , 'LineWidth' , 2 );
          drawnow;
       end
    end

    % stopping criteria - - - - - - - - - - - - - - - - - - - - - - - - - -

    if ng <= eps * ng0
       status = 'optimal';
       fprintf( '\n' );
       break;
    end
    
    if feval > MaxFeval
       status = 'stopped';
       fprintf( '\n' );
       break;
    end

    % compute search direction- - - - - - - - - - - - - - - - - - - - - - -
    % formulae could be streamlined somewhat and some norms could be saved
    % from previous iterations

    if sgiter     % iteration without previous d == standard gradient
       d = - lastg;     % now the previous d is defined
       sgiter = false;  % next iteration will be a "normal" one
       fprintf( '\t********' );
    else          % normal iteration, use appropriate NCG formula
       if rstart > 0 && mod( iter , n * rstart ) == 0
          beta = 0;  % ... unless a restart is being performed
          fprintf( '\t(res)   ' );
       else
          switch wf
             case 0     % Fletcher-Reeves
                beta = ( ng / norm( pastg ) )^2;
             case 1     % Polak-Ribiere
                beta = ( lastg' * ( lastg - pastg ) ) / norm( pastg )^2;
                beta = max( [ beta 0 ] );
             case 2     % Hestenes-Stiefel
                beta = ( lastg' * ( lastg - pastg ) ) / ...
                       ( ( lastg - pastg )' * pastd );
                beta = max( [ beta 0 ] );
             otherwise  % Dai-Yuan
                beta = ng^2 / ( ( lastg - pastg )' * pastd );
          end
          fprintf( '\t%1.2e' , beta );
       end

       if beta ~= 0
          d = - lastg + beta * pastd;
       else
          d = - lastg;
       end
    end

    pastg = lastg;  % previous gradient
    pastd = d;      % previous search direction
    
    % compute step size - - - - - - - - - - - - - - - - - - - - - - - - - -

    phip0 = lastg' * d;

    if AWLS
       [ a , v ] = ArmijoWolfeLS( v , phip0 , astart , m1 , m2 , tau );
    else
       [ a , v ] = BacktrackingLS( v , phip0 , astart , m1 , tau );
    end

    % output statistics - - - - - - - - - - - - - - - - - - - - - - - - - -

    fprintf( '\t   %1.2e\n' , a );
 
    if a <= abs( mina )
       if mina >= 0 || beta == 0
          % if beta was 0, then the last direction used was already the
          % "pure" anti-gradient, hence failure in the line search cannot
          % be attributed to a wrong value of beta, and even more
          % importantly it cannot be solved by a reset
          status = 'error';
          fprintf( '\n' );
          break;
       else
          % in case of problems with the direction, force a restart
          sgiter = true;
       end
    end
    
    if v <= MInf
       status = 'unbounded';
       fprintf( '\n' );
       break;
    end

    % compute new point - - - - - - - - - - - - - - - - - - - - - - - - - -

    % possibly plot the trajectory
    if n == 2 && Plotf == 1
       PXY = [ x ,  lastx ];
       line( 'XData' , PXY( 1 , : ) , 'YData' , PXY( 2 , : ) , ...
             'LineStyle' , '-' , 'LineWidth' , 2 ,  'Marker' , 'o' , ...
             'Color' , [ 0 0 0 ] );
       drawnow;
    end
    
    x = lastx;
    ng = norm( lastg );

    % iterate - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    iter = iter +1;

    if Interactive
       pause;
    end
end

% end of main loop- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% inner functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function [ phi , varargout ] = f2phi( alpha )
%
% computes and returns the value of the tomography at alpha
%
%    phi( alpha ) = f( x + alpha * d )
%
% if Plotf > 2 saves the data in gap() for plotting
%
% if the second output parameter is required, put there the derivative
% of the tomography in alpha
%
%    phi'( alpha ) = < \nabla f( x + alpha * d ) , d >
%
% saves the point in lastx, the gradient in lastg and increases feval

   lastx = x + alpha * d;
   [ phi , lastg ] = f( lastx );

   if Plotf > 2
      if fStar > - Inf
         gap( end + 1 ) = ( phi - fStar ) / max( [ abs( fStar ) 1 ] );
      else
         gap( end + 1 ) = phi;
      end
   end
   
   if nargout > 1
      varargout{ 1 } = d' * lastg;
   end
   
   feval = feval + 1;
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function [ a , phia ] = ArmijoWolfeLS( phi0 , phip0 , as , m1 , m2 , tau )

% performs an Armijo-Wolfe Line Search.
%
% phi0 = phi( 0 ), phip0 = phi'( 0 ) < 0
%
% as > 0 is the first value to be tested: if phi'( as ) < 0 then as is
% divided by tau < 1 (hence it is increased) until this does not happen
% any longer
%
% m1 and m2 are the standard Armijo-Wolfe parameters; note that the strong
% Wolfe condition is used
%
% returns the optimal step and the optimal f-value

lsiter = 1;  % count iterations of first phase
while feval <= MaxFeval 
   [ phia , phips ] = f2phi( as );

   if ( phia <= phi0 + m1 * as * phip0 ) && ...
      ( abs( phips ) <= - m2 * phip0 )
      fprintf( ' %2d   ' , lsiter );
      a = as;
      return;  % Armijo + strong Wolfe satisfied, we are done

   end
   if phips >= 0  % derivative is positive, break
      break;
   end
   as = as / tau;
   lsiter = lsiter + 1;
end    

fprintf( ' %2d ' , lsiter );
lsiter = 1;  % count iterations of second phase

am = 0;
a = as;
phipm = phip0;
while ( feval <= MaxFeval ) && ( ( as - am ) ) > abs( mina ) && ...
      ( phips > 1e-12 )

   % compute the new value by safeguarded quadratic interpolation
   a = ( am * phips - as * phipm ) / ( phips - phipm );
   a = max( [ am + ( as - am ) * sfgrd ...
            min( [ as - ( as - am ) * sfgrd  a ] ) ] );

   % compute phi( a )
   [ phia , phip ] = f2phi( a );

   if ( phia <= phi0 + m1 * a * phip0 ) && ( abs( phip ) <= - m2 * phip0 )
      break;  % Armijo + strong Wolfe satisfied, we are done
   end

   % restrict the interval based on sign of the derivative in a
   if phip < 0
      am = a;
      phipm = phip;
   else
      as = a;
      if as <= abs( mina )
         break;
      end
      phips = phip;
   end
   lsiter = lsiter + 1;
end

fprintf( '%2d' , lsiter );

end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function [ as , phia ] = BacktrackingLS( phi0 , phip0 , as , m1 , tau )

% performs a Backtracking Line Search.
%
% phi0 = phi( 0 ), phip0 = phi'( 0 ) < 0
%
% as > 0 is the first value to be tested, which is decreased by
% multiplying it by tau < 1 until the Armijo condition with parameter
% m1 is satisfied
%
% returns the optimal step and the optimal f-value

lsiter = 1;  % count ls iterations
while feval <= MaxFeval && as > abs( mina )
   [ phia , ~ ] = f2phi( as );
   if phia <= phi0 + m1 * as * phip0  % Armijo satisfied
      break;                          % we are done
   end
   as = as * tau;
   lsiter = lsiter + 1;
end

fprintf( '\t%2d' , lsiter );

end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end  % the end- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




