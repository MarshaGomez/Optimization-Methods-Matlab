function [ x , status ] =  HBG( f , varargin )

%function [ x , status ] = HBG( f , x , alpha , beta , eps , MaxIter ,
%                               MInf )
%
% Apply a Heavy Ball Gradient approach for the minimization of the
% provided function f, which must have the following interface:
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
% - alpha (real scalar, optional, default value 1): the fixed stepsize of
%   the Heavy Ball Gradient approach (along the anti-gradient).
%
% - beta (real scalar, optional, default value 0.9): the fixed weight of
%   the momentum term
%
%        beta * || x^i - x^{i - 1} ||
%
%   Note that beta has to be >= 0, although 0 is accepted which turns the
%   Heavy Ball Gradient approach into a "Light" Ball Gradient approach,
%   i.e., a standard Gradient approach with fixed stepsize.
%
% - eps (real scalar, optional, default value 1e-6): the accuracy in the 
%   stopping criterion: the algorithm is stopped when the norm of the
%   gradient is less than or equal to eps. If a negative value is provided,
%   this is used in a *relative* stopping criterion: the algorithm is
%   stopped when the norm of the gradient is less than or equal to
%   (- eps) * || norm of the first gradient ||.
%
% - MaxIter (integer scalar, optional, default value 300): the maximum
%   number of iterations == function evaluations.
%
% - MInf (real scalar, optional, default value -Inf): if the algorithm
%   determines a value for f() <= MInf this is taken as an indication that
%   the problem is unbounded below and computation is stopped
%   (a "finite -Inf").
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
 Date: 10-11-22
 Version 1.01
 Copyright Antonio Frangioni
 =======================================
%}

Plotf = 1;
% 0 = nothing is plotted
% 1 = the level sets of f and the trajectory are plotted (when n = 2)
% 2 = the function value / gap are plotted

if Plotf == 2
   gap = [];
end

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
   alpha = varargin{ 2 };
   if ~ isscalar( alpha )
      error( 'beta is not a real scalar' );
   end
   if alpha <= 0
      error( 'alpha must be positive' );
   end
else
   alpha = 1;
end

if length( varargin ) > 2
   beta = varargin{ 3 };
   if ~ isscalar( beta )
      error( 'beta is not a real scalar' );
   end
   if beta < 0
      error( 'beta must be non-negative' );
   end
else
   beta = 0.9;
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
   MaxIter = round( varargin{ 5 } );
   if ~ isscalar( MaxIter )
      error( 'MaxFeval is not an integer scalar' );
   end
else
   MaxIter = 1000;
end

if length( varargin ) > 5
   MInf = varargin{ 6 };
   if ~ isscalar( MInf )
      error( 'MInf is not a real scalar' );
   end
else
   MInf = - Inf;
end

% initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

fprintf( 'Heavy Ball Gradient method\n');
if fStar > - Inf
   fprintf( 'feval\trel gap\t\tbest gap');
else
   fprintf( 'feval\tf(x)\tfbest');
end
fprintf( '\t|| g(x) ||\n\n');

if Plotf == 2
   xlim( [ 0 MaxIter ] );
   ax = gca;
   ax.FontSize = 16;
   ax.Position = [ 0.03 0.07 0.95 0.92 ];
   ax.Toolbar.Visible = 'off';
end

[ v , g ] = f( x );
ng = norm( g );
vbest = v;
if eps < 0
   ng0 = - ng;  % norm of first subgradient: why is there a "-"? ;-)
else
   ng0 = 1;     % un-scaled stopping criterion
end

pastd = zeros( n , 1 );  % the direction at the previous iteration
feval = 1;               % f() evaluations count

% main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

while true
 
    % output statistics - - - - - - - - - - - - - - - - - - - - - - - - - -
  
    if fStar > - Inf
       gapk = ( v - fStar ) / max( [ abs( fStar ) 1 ] );
       bstgapk = ( vbest - fStar ) / max( [ abs( fStar ) 1 ] );
       
       fprintf( '%4d\t%1.4e\t%1.4e\t%1.4e\n' , feval , gapk , ...
                bstgapk , ng );

       if Plotf == 2
          gap( end + 1 ) = gapk;
          semilogy( gap , 'Color' , 'k' , 'LineWidth' , 2 );
          ylim( [ 1e-15 1e+1 ] );
          drawnow;
       end
    else
       fprintf( '%4d\t%1.8e\t%1.8e\t\t%1.4e\n' , feval , v , vbest , ng );

       if Plotf == 2
          gap( end + 1 ) = v;
          plot( gap , 'Color' , 'k' , 'LineWidth' , 2 );
          drawnow;
       end
    end
    
    % stopping criteria - - - - - - - - - - - - - - - - - - - - - - - - - -

    if ng <= eps * ng0
       status = 'optimal';
       break;
    end
    
    if feval > MaxIter
       status = 'stopped';
       break;
    end
    
    if v <= MInf
       status = 'unbounded';
       break;
    end

    % compute deflected gradient direction- - - - - - - - - - - - - - - - -
    
    d = - alpha * g + beta * pastd;     
 
    % compute new point - - - - - - - - - - - - - - - - - - - - - - - - - -

    % possibly plot the trajectory
    if n == 2 && Plotf == 1
       PXY = [ x ,  x + d ];
       line( 'XData' , PXY( 1 , : ) , 'YData' , PXY( 2 , : ) , ...
             'LineStyle' , '-' , 'LineWidth' , 2 ,  'Marker' , 'o' , ...
             'Color' , [ 0 0 0 ] );
       drawnow;
    end

    x = x + d;
    pastd = d;
  
    % compute f() - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    [ v , g ] = f( x );
    ng = norm( g );
    if v < vbest
       vbest = v;
    end
    feval = feval + 1;

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

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end  % the end- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




