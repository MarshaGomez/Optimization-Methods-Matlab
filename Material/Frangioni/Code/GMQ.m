function [ x , status ] =  GMQ( Q , q , varargin )

%function [ x , status ] = GMQ( Q , q , x , fStar , alpha , MaxIter , eps )
%
% Apply the Gradient Method (a.k.a., Steepest Descent algorithm) to the
% minimization of the quadratic function
%
%   f( x ) = 1/2 x^T Q x + q x
%
% Input:
%
% - Q ([ n x n ] real symmetric matrix, not necessarily positive
%   semidefinite): the quadratic part of f
%
% - q ([ n x 1 ] real column vector): the linear part of f
%
% - x ([ n x 1 ] real column vector or empty, optional): the point where to
%   start the algorithm from; if not provided or empty, the all-0 n-vector
%   is used
%
% - fStar (real scalar, optional, default value Inf): optimal value of f.
%   if a non-Inf value is provided it is used to print out stasistics about
%   the convergence speed of the algorithm
%
% - alpha (real scalar, optional, default value 0): if alpha > 0, then the
%   fixed-stepsize version of the algorithm is run with alpha as the fixed
%   stepsize, otherwise the standard exact line search is used
%
% - MaxIter (integer scalar, optional, default value 1000): the maximum
%   number of iterations
%
% - eps (real scalar, optional, default value 1e-6): the accuracy in the 
%   stopping criterion: the algorithm is stopped when the norm of the
%   gradient is less than or equal to eps
%
% Output:
%
% - x ([ n x 1 ] real column vector): either the best solution found so far
%   (possibly the optimal one) or a direction proving the problem is
%   unbounded below, depending on case
%
% - status (string): a string describing the status of the algorithm at
%   termination 
%
%   = 'optimal': the algorithm terminated having proven that x is a(n
%     approximately) optimal solution, i.e., the norm of the gradient at x
%     is less than the required threshold
%
%   = 'unbounded': the algorithm terminated having proven that the problem
%     is unbounded below: x contains a direction along which f is
%     decreasing to - Inf, either because f is linear along x and the
%     directional derivative is not zero, or because x is a direction with
%     negative curvature
%
%   = 'stopped': the algorithm terminated having exhausted the maximum
%     number of iterations: x is the best solution found so far, but not
%     necessarily the optimal one
%
%{
 =======================================
 Author: Antonio Frangioni
 Date: 26-12-22
 Version 0.2
 Copyright Antonio Frangioni
 =======================================
%}

Plotf = 0;
% 0 = nothing is plotted
% 1 = the function value / gap are plotted
% 2 = the level sets of f and the trajectory are plotted (when n = 2)

Interactive = false;  % if we pause at every iteration

Streamlined = true;  % if the streamlined version of the algorithm, with
                     % only one O( n^2 ) operation per iteration, is used

% reading and checking input- - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if ~ isreal( Q )
   error( 'Q not a real matrix' );
end

n = size( Q , 1 );

if n <= 1
   error( 'Q is too small' );
end

if n ~= size( Q , 2 )
   error( 'Q is not square' );
end

if ~ isreal( q )
   error( 'q not a real vector' );
end

if size( q , 2 ) ~= 1
   error( 'q is not a (column) vector' );
end

if size( q , 1 ) ~= n
   error( 'q size does not match with Q' );
end

if isempty( varargin ) || isempty( varargin{ 1 } )
   x = zeros( n , 1 );
else
   x = varargin{ 1 };

   if ~ isreal( x )
      error( 'x not a real vector' );
   end

   if size( x , 2 ) ~= 1
      error( 'x is not a (column) vector' );
   end

   if size( x , 1 ) ~= n
      error( 'x size does not match with Q' );
   end
end

if length( varargin ) > 1
   fStar = varargin{ 2 };
   if ~ isreal( fStar ) || ~ isscalar( fStar )
      error( 'fStar is not a real scalar' );
   end
else
   fStar = - Inf;
end

if length( varargin ) > 2
   alpha = varargin{ 3 };
   if ~ isreal( alpha ) || ~ isscalar( alpha )
      error( 'alpha is not a real scalar' );
   end
else
   alpha = 0;
end

if length( varargin ) > 3
   MaxIter = round( varargin{ 4 } );
   if ~ isscalar( MaxIter )
      error( 'MaxIter is not an integer scalar' );
   end
   if MaxIter < 1
      error( 'MaxIter too small' );
   end
else
   MaxIter = 1000;
end

if length( varargin ) > 4
   eps = varargin{ 5 };
   if ~ isreal( eps ) || ~ isscalar( eps )
      error( 'eps is not a real scalar' );
   end
   if eps < 0
      error( 'eps can not be negative' );
   end
else
   eps = 1e-6;
end

% initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

fprintf( 'Gradient method for quadratic functions ' );
if alpha == 0
   fprintf( '(optimal stepsize)\n' );
else
   fprintf( '(fixed stepsize)\n' );
end
fprintf( 'iter\tf(x)\t\t\t||g||');
if fStar > - Inf
   fprintf( '\t\tgap\t\trate');
   prevf = Inf;
end
if alpha == 0
   fprintf( '\t\talpha' );
end
fprintf( '\n\n' );

i = 0;

if Plotf == 1
   gap = [];          
end

if Streamlined
   g = Q * x + q;
end

% main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

while true
    % compute function value and gradient - - - - - - - - - - - - - - - - -

    if ~ Streamlined
       g = Q * x + q; 
    end
    
    ng = norm( g );
    f = ( g + q )' * x / 2;  % 1/2 x^T Q x + q x = 1/2 ( x^T Q x + 2 q x )
                             % = 1/2 x^T ( Q x + q + q ) = 1/2 ( q + g ) x
    i = i + 1;

    % output statistics - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
    fprintf( '%4d\t%1.8e\t\t%1.4e' , i , f , ng );
    if fStar > - Inf
       gapk = ( f - fStar ) / max( [ abs( fStar ) 1 ] );

       fprintf( '\t%1.4e' , gapk );
       if prevf < Inf
          fprintf( '\t%1.4e' , ( f - fStar ) / ( prevf - fStar ) );
       else
          fprintf( '\t\t' );
       end
       prevf = f;

       if Plotf == 1
          gap( end + 1 ) = gapk;
          semilogy( gap , 'Color' , 'k' , 'LineWidth' , 2 );
          xlim( [ 0 MaxIter ] );
          ylim( [ 1e-15 inf ] );
          ax = gca;
          ax.FontSize = 16;
          ax.Position = [ 0.03 0.07 0.95 0.92 ];
          ax.Toolbar.Visible = 'off';
       end
    end
   
    % stopping criteria - - - - - - - - - - - - - - - - - - - - - - - - - -
    if ng <= eps
       status = 'optimal';
       if alpha == 0
          fprintf( '\n' );
       end
       break;
    end
    
    if i > MaxIter
       status = 'stopped';
       if alpha == 0
          fprintf( '\n' );
       end
       break;
    end
    
    % compute step size - - - - - - - - - - - - - - - - - - - - - - - - - -
    % meanwhile, check if f is unbounded below
    % note that if alpha > 0 this is only used for the unboundedness check
    % which is a bit of a waste, but there you go; anyway, in the
    % streamlined version this only costs O( n )

    if Streamlined
       v = Q * g; 
       den = g' * v;
    else
       den = g' * Q * g;
    end

    if den <= 1e-14
       % this is actually two different cases:
       %
       % - g' * Q * g = 0, i.e., f is linear along g, and since the
       %   gradient is not zero, it is unbounded below
       %
       % - g' * Q * g < 0, i.e., g is a direction of negative curvature for
       %   f, which is then necessarily unbounded below
       %
       if alpha == 0
          fprintf( '\n' );
       end
       fprintf( 'g'' * Q * g = %1.4e ==> unbounded\n' , den );
       status = 'unbounded';
       break;
    end
    
    if alpha > 0
       t = alpha;
    else
       t = ng^2 / den;  % stepsize
       fprintf( '\t%1.2e' , t );
    end

    fprintf( '\n' );

    % compute new point - - - - - - - - - - - - - - - - - - - - - - - - - -

    % possibly plot the trajectory
    if n == 2 && Plotf == 2
       PXY = [ x ,  x - t * g ];
       line( 'XData' , PXY( 1 , : ) , 'YData' , PXY( 2 , : ) , ...
             'LineStyle' , '-' , 'LineWidth' , 2 ,  'Marker' , 'o' , ...
             'Color' , [ 0 0 0 ] );
    end
    
    x = x - t * g;

    if Streamlined
       g = g - t * v;
    end
    
    % iterate - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if Interactive
       pause;
    end
end

% end of main loop- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end  % the end- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




