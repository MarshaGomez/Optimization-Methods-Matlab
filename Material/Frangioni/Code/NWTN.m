function [ x , status ] = NWTN( f , varargin )

%function [ x , status ] = NWTN( f , x , eps , MaxFeval , m1 , m2 , delta ,
%                                tau , sfgrd , MInf , mina )
%
% Apply a classical Newton's method for the minimization of the provided
% function f, which must have the following interface:
%
%   [ v , g , H ] = f( x )
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
% - H (real, [ n x n ] real matrix) must only be specified if x ~= [],
%   and it is the Hessian of f() at x. If no such information is
%   available, the function throws error.
%
% The other [optional] input parameters are:
%
% - x (either [ n x 1 ] real vector or [], default []): starting point.
%   If x == [], the default starting point provided by f() is used.
%
% - eps (real scalar, optional, default value 1e-6): the accuracy in the 
%   stopping criterion: the algorithm is stopped when the norm of the
%   gradient is less than or equal to eps. If a negative value is provided,
%   this is used in a *relative* stopping criterion: the algorithm is
%   stopped when the norm of the gradient is less than or equal to
%   (- eps) * || norm of the first gradient ||.
%
% - MaxFeval (integer scalar, optional, default value 1000): the maximum
%   number of function evaluations (hence, iterations will be not more than
%   MaxFeval because at each iteration at least a function evaluation is
%   performed, possibly more due to the line search).
%
% - m1 (real scalar, optional, default value 1e-4): first parameter of the
%   Armijo-Wolfe-type line search (sufficient decrease). Has to be in (0,1)
%
% - m2 (real scalar, optional, default value 0.9): typically the second
%   parameter of the Armijo-Wolfe-type line search (strong curvature
%   condition). It should to be in (0,1); if not, it is taken to mean that
%   the simpler Backtracking line search should be used instead
%
% - delta (real scalar, optional, default value 1e-6): minimum positive
%   value for the eigenvalues of the modified Hessian used to compute the
%   Newton direction
%
% - tau (real scalar, optional, default value 0.9): scaling parameter for
%   the line search. In the Armijo-Wolfe line search it is used in the
%   first phase: if the derivative is not positive, then the step is
%   divided by tau (which is < 1, hence it is increased). In the 
%   Backtracking line search, each time the step is multiplied by tau
%   (hence it is decreased).
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
% - MInf (real scalar, optional, default value -Inf): if the algorithm
%   determines a value for f() <= MInf this is taken as an indication that
%   the problem is unbounded below and computation is stopped
%   (a "finite -Inf").
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
 Date: 29-10-21
 Version 1.21
 Copyright Antonio Frangioni
 =======================================
%}

Plotf = 1;
% 0 = nothing is plotted
% 1 = the level sets of f and the trajectory are plotted (when n = 2)
% 2 = the function value / gap are plotted, iteration-wise
% 3 = the function value / gap are plotted, function-evaluation-wise

Interactive = true;  % if we pause at every iteration

if Plotf > 1
   if Plotf == 2
      MaxIter = 50;  % expected number of iterations for the gap plot
   else
      MaxIter = 70;  % expected number of iterations for the gap plot
   end    
   gap = [];

   xlim( [ 0 MaxIter ] );
   ax = gca;
   ax.FontSize = 16;
   ax.Position = [ 0.03 0.07 0.95 0.92 ];
   ax.Toolbar.Visible = 'off';
end

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
   m1 = 1e-4;
end

if length( varargin ) > 4
   m2 = varargin{ 5 };
   if ~ isscalar( m1 )
      error( 'm2 is not a real scalar' );
   end
else
   m2 = 0.9;
end

AWLS = ( m2 > 0 && m2 < 1 );

if length( varargin ) > 5
   delta = varargin{ 6 };
   if ~ isscalar( delta )
      error( 'delta is not a real scalar' );
   end
   if delta < 0
      error( 'delta must be > 0' );
   end       
else
   delta = 1e-6;
end

if length( varargin ) > 6
   tau = varargin{ 7 };
   if ~ isscalar( tau )
      error( 'tau is not a real scalar' );
   end
   if tau <= 0 || tau >= 1
      error( 'tau is not in (0 ,1)' );
   end       
else
   tau = 0.9;
end

if length( varargin ) > 7
   sfgrd = varargin{ 8 };
   if ~ isscalar( sfgrd )
      error( 'sfgrd is not a real scalar' );
   end
   if sfgrd <= 0 || sfgrd >= 1
      error( 'sfgrd is not in (0, 1)' );
   end
else
   sfgrd = 0.2;
end

if length( varargin ) > 8
   MInf = varargin{ 9 };
   if ~ isscalar( MInf )
      error( 'MInf is not a real scalar' );
   end
else
   MInf = - Inf;
end

if length( varargin ) > 9
   mina = varargin{ 10 };
   if ~ isscalar( mina )
      error( 'mina is not a real scalar' );
   end
   if mina < 0
      error( 'mina is < 0' );
   end       
else
   mina = 1e-16;
end

% "global" variables- - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

lastx = zeros( n , 1 );  % last point visited in the line search
lastg = zeros( n , 1 );  % gradient of lastx
lastH = zeros( n , n );  % Hessian of lastx
d = zeros( n , 1 );      % Newton's direction
feval = 0;               % f() evaluations count ("common" with LSs)

% initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

fprintf( 'Newton''s method\n');
if fStar > - Inf
   fprintf( 'feval\trel gap\t\t|| g(x) ||\trate\t\tdelta');
   prevv = Inf;
else
   fprintf( 'feval\tf(x)\t\t\t|| g(x) ||\tdelta');
end
fprintf( '\tls it\ta*');
fprintf( '\n\n' );

v = f2phi( 0 );
ng = norm( lastg );
if eps < 0
   ng0 = - ng;  % norm of first subgradient: why is there a "-"? ;-)
else
   ng0 = 1;     % un-scaled stopping criterion
end

% main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

while true

    % output statistics - - - - - - - - - - - - - - - - - - - - - - - - - -
  
    if fStar > - Inf
       gapk = ( v - fStar ) / max( [ abs( fStar ) 1 ] );

       fprintf( '%4d\t%1.4e\t%1.4e' , feval , gapk , ng );
       if prevv < Inf
          fprintf( '\t%1.4e' , ( v - fStar ) / ( prevv - fStar ) );
       else
          fprintf( '\t\t' );           
       end
       prevv = v;

       if Plotf > 1
          if Plotf >= 2
             gap( end + 1 ) = gapk;
          end
          semilogy( gap , 'Color' , 'k' , 'LineWidth' , 2 );
          if Plotf == 2
             ylim( [ 1e-15 1e+1 ] );
          else
             ylim( [ 1e-15 1e+4 ] );
          end
          drawnow;
       end
    else
       fprintf( '%4d\t%1.8e\t\t%1.4e' , feval , v , ng );

       if Plotf > 1
          if Plotf >= 2
             gap( end + 1 ) = v;
          end
          plot( gap , 'Color' , 'k' , 'LineWidth' , 2 );
          drawnow;
       end
    end
    
    % stopping criteria - - - - - - - - - - - - - - - - - - - - - - - - - -

    if ng <= eps * ng0
       status = 'optimal';
       break;
    end
    
    if feval > MaxFeval
       status = 'stopped';
       break;
    end

    % compute Newton's direction- - - - - - - - - - - - - - - - - - - - - -

    lambdan = eigs( lastH , 1 , 'sa' );  % smallest eigenvalue
    if lambdan < delta
       fprintf( '\t%1.2e' , delta - lambdan );
       lastH = lastH + ( delta - lambdan ) * eye( n );
    else
       fprintf( '\t0.00e+00' );
    end
    d = - lastH \ lastg;
    phip0 = lastg' * d;

    % compute step size - - - - - - - - - - - - - - - - - - - - - - - - - -
    % in Newton's method, the default initial stepsize is 1
    
    if AWLS
       [ a , v ] = ArmijoWolfeLS( v , phip0 , 1 , m1 , m2 , tau );
    else
       [ a , v ] = BacktrackingLS( v , phip0 , 1 , m1 , tau );
    end

    % output statistics - - - - - - - - - - - - - - - - - - - - - - - - - -

    fprintf( '\t%1.2e' , a );
    fprintf( '\n' );

    if a <= mina
       status = 'error';
       break;
    end
    
    if v <= MInf
       status = 'unbounded';
       break;
    end

    % compute new point - - - - - - - - - - - - - - - - - - - - - - - - - -

    % possibly plot the trajectory
    if n == 2 && Plotf == 1
       PXY = [ x , lastx ];
       line( 'XData' , PXY( 1 , : ) , 'YData' , PXY( 2 , : ) , ...
             'LineStyle' , '-' , 'LineWidth' , 2 ,  'Marker' , 'o' , ...
             'Color' , [ 0 0 0 ] );
    end
    
    x = lastx;
    ng = norm( lastg );

    % iterate - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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
% saves the point in lastx, the gradient in lastg, the Hessian in lasth,
% and increases feval

   lastx = x + alpha * d;
   [ phi , lastg , lastH ] = f( lastx );

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

   if ( phia <= phi0 + m1 * as * phip0 ) && ( abs( phips ) <= - m2 * phip0 )
      fprintf( '  %2d' , lsiter );
      a = as;
      return;  % Armijo + strong Wolfe satisfied, we are done

   end
   if phips >= 0
      break;
   end
   as = as / tau;
   lsiter = lsiter + 1;
end    

fprintf( '  %2d ' , lsiter );
lsiter = 1;  % count iterations of second phase

am = 0;
a = as;
phipm = phip0;
while ( feval <= MaxFeval ) && ( ( as - am ) ) > mina && ( phips > 1e-12 )

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
      if as <= mina
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
while feval <= MaxFeval && as > mina
   [ phia , ~ ] = f2phi( as );
   if phia <= phi0 + m1 * as * phip0  % Armijo satisfied
      break;                        % we are done
   end
   as = as * tau;
   lsiter = lsiter + 1;
end

fprintf( '  %2d' , lsiter );

end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end  % the end- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




