function [ x , status ] =  DIS( f , varargin )

%function [ x , status ] = DIS( f , range , sfgrd , delta , MaxFeval )
%
% Apply the classical Dichotomic Search for the minimization of the
% provided one-dimensional function f, which must have the following
% interface:
%
%   [ v , varargout ] = f( x )
%
% Input:
%
% - x is either a scalar real denoting the input of f(), or [] (empty).
%
% Output:
%
% - v (real, scalar): if x == [] this is the best known lower bound on
%   the global optimum of f() on the standard interval in which f() is
%   supposed to be minimised (see next). If x ~= [] then v = f(x).
%
% - g (real, either scalar or a [ 2 x 1 ] matrix denoting an interval) is
%   the first optional argument. This also depends on x. if x == [] then
%   g is a [ 2 x 1 ] matrix denoting the standard interval in which f()
%   is supposed to be minimised (into which v is the minimum). f() is
%   never called with x ~= [].
%
% IMPORTANT NOTE: the function requires f() to be able to provide the first
%                 derivative, and it requires that the interval
%  [ x_- , x_+ ] is chosen such that f'( x_- ) <= 0 and f'( x_+ ) >= 0.
%
% The Dichotomic Search can either be "blind" (new point right in the
% middle of the interval) or using a safeguarded quadratic Interpolation
% to choose it, as dictated by the other [optional] input parameters:
%
% - range: (either [ 2 x 1 ] real vector or [], default []): the range
%   in which the local minimum has to be seached; if range == [], the
%   default range point provided by f() is used.
%
% - sfgrd (real scalar, default value 0): if sfgrd == 0, the Dichotomic
%   Search is "blind", i .e., the new point is always chosen right in the
%   middle of the current interval. Otherwise, it must be 0 < sfgrd < 0.5
%   and a safeguarded quadratic Interpolation technique is used to choose
%   it, where it is guaranteed that at least ( 1 - sfgrd ) of the current
%   interval will be discarded.
%
% - eps (real scalar, default value 1e-6): the accuracy in the stopping
%   criterion: the algorithm is stopped when a point is found such that
%   the absolute value of the derivative is less than or equal to eps.
%
% - MaxFeval (integer scalar, default value 100): the maximum number of
%   function evaluations (hence, iterations will be not more than
%   MaxFeval - 2 because at each iteration one function evaluation is
%   performed, except in the first one when two are).
%
% Output:
%
% - x (real scalar): the best solution found so far.
%
% - status (string): a string describing the status of the algorithm at
%   termination 
%
%   = 'optimal': the algorithm terminated having proven that x is a(n
%     approximately) optimal solution, i.e., the diameter of the
%     restricted range is less than or equal to delta.
%
%   = 'empty': the provided range is empty (x_- > x_+) and therefore
%     such is the optimization problem
%
%   = 'stopped': the algorithm terminated having exhausted the maximum
%     number of iterations: x is the best solution found so far, but not
%     necessarily the optimal one
%
%   = 'error': the algorithm found a numerical error that prevents it from
%     continuing optimization
%
% TODO: implement a warm-op phase whereby if f'( x_- ) > 0 then x_- is
%       "quickly moved left" until the derivative is negative, and,
%       symmetrically, if f'( x_+ ) < 0 then it is "quickly moved right".
%
%{
 =======================================
 Author: Antonio Frangioni
 Date: 09-08-21
 Version 0.10
 Copyright Antonio Frangioni
 =======================================
%}

Plotg = 1;
% 0 = nothing is plotted
% 1 = the function value / gap are plotted
% 2 = the function and the model (if used) are plotted

Interactive = true;  % if we pause at every iteration

% reading and checking input- - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if ~ isa( f , 'function_handle' )
   error( 'f not a function' );
end

if isempty( varargin ) || isempty( varargin{ 1 } )
   [ fStar , range ] = f( [] );
else
   fStar = - Inf;  % if the range is not the standard one, we can't trust
                   % the standard global minima
   range = varargin{ 1 };
end

if ~ isreal( range )
   error( 'range not a real vector' );
end

if ( size( range , 1 ) ~= 1 ) || ( size( range , 2 ) ~= 2 )
   error( 'range is not a [ 1 x 2 ] vector' );
end
   
xm = range( 1 );  % x_-
xp = range( 2 );  % x_+

if xm > xp
   fprintf( 'range is empty\n' );
   x = NaN;
   status = 'empty';
   return;
end

[ fxm , f1xm ] = f( xm );
if( f1xm > 0 )
    error( 'f''( x_- ) must be <= 0' );
end

[ fxp , f1xp ] = f( xp );
if( f1xp < 0 )
    error( 'f''( x_+ ) must be >= 0' );
end

if length( varargin ) > 1
   sfgrd = varargin{ 2 };
   if ~ isreal( sfgrd ) || ~ isscalar( sfgrd )
      error( 'sfgrd is not a real scalar' );
   end
   if ( sfgrd < 0 ) || ( sfgrd >= 0.5 )
      error( 'sfgrd must be in [ 0 , 1/2 )' );
   end
else
   sfgrd = 0;
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
   MaxFeval = round( varargin{ 4 } );
   if ~ isscalar( MaxFeval )
      error( 'MaxFeval is not an integer scalar' );
   end
else
   MaxFeval = 100;
end

% initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

feval = 2;

status = 'optimal';

if f1xm >= - eps
   x = xm;
   return;
end

if f1xp <= eps
   x = xp;
   return;
end

fbest = min( [ fxp fxm ] );
f1x = min( [ -f1xp f1xm ] );

if Plotg == 1
   gap = [];          
end

if sfgrd == 0
   fprintf( 'Dichotomic search\n');
else
   fprintf( ...
    'Dichotomic search with safeguarded interpolation (%1.4f)\n' , sfgrd );
end
if fStar > - Inf
   fprintf( 'feval\trel gap\t\tx_-\t\tx_+\t\tx\t\tf''(x)\n');
else
   fprintf( 'feval\tfbest\t\tx_-\t\tx_+\t\tx\t\tf''(x)\n');
end

% main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

while true

    % stopping criteria - - - - - - - - - - - - - - - - - - - - - - - - - -

    if feval > MaxFeval
       status = 'stopped';
       break;
    end
    
    % main logic- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if sfgrd == 0
       % compute the new point as the middle of the inteval
       x = ( xm + xp ) / 2;
    else
       % compute the new point by safeguarded quadratic interpolation
       sfxp = xp - ( xp - xm ) * sfgrd;
       sfxm = xm + ( xp - xm ) * sfgrd;
       x = ( xm * f1xp - xp * f1xm ) / ( f1xp - f1xm );
       x = max( [ sfxm min( [ sfxp x ] ) ] );
    end

    [ fx , f1x ] = f( x );  % compute f( x ) and f'( x )
 
    if fx < fbest
       fbest = fx;
    end
    
    feval = feval + 1;

    % output statistics - - - - - - - - - - - - - - - - - - - - - - - - - -
  
    if fStar > - Inf
       gapk = ( fbest - fStar ) / max( [ abs( fStar ) 1 ] );

       if Plotg == 1
          gap( end + 1 ) = gapk;
          semilogy( gap , 'Color' , 'k' , 'LineWidth' , 2 );
          xlim( [ 0 35 ] );
          ylim( [ 1e-15 inf ] );
          ax = gca;
          ax.FontSize = 16;
          ax.Position = [ 0.03 0.07 0.95 0.92 ];
          ax.Toolbar.Visible = 'off';
       end
     else
       gapk = fbest;
    end
    fprintf( '%4d\t%1.4e\t%1.8e\t%1.8e\t%1.8e\t%1.4e\n' , feval , ...
             gapk , xm , xp , x , f1x );

    if Plotg == 2
       xmp = xm - ( xp - xm ) / 20;
       xpp = xp + ( xp - xm ) / 20;

       warning( 'off' , 'all' );
       fplot( @(x) f( x ) , [ xmp xpp ] , 'Color' , 'k' , ...
              'LineWidth' , 1 );

       xlim( [ xmp xpp ] );
       yticks( [] );
       ax = gca;
       ax.FontSize = 16;
       ax.Toolbar.Visible = 'off';
       ax.Position = [ 0.025 0.05 0.95 0.95 ];
       if sfgrd ~= 0
          hold on;
          a = ( f1xp - f1xm ) / ( 2 * ( xp - xm ) );
          b = ( xp * f1xm - xm * f1xp ) / ( xp - xm );
          % a xm^2 + b xm + c == fxm   ==>
          % c == fxm - a xm^2 - b xm
          c = fxm - a * xm^2 - b * xm;
          fplot( @(x) a * x^2 + b * x + c , [ xmp xpp ] , ... 
                'Color' , 'b' , 'LineWidth' , 1 );
          xticks( [ xmp xm sfxm sfxp xp xpp ] );
          xticklabels( { num2str( xmp , '%1.1g' ) , 'x-' , 'sx-' , ...
                         'sx+'  , 'x+' , num2str( xpp , '%1.1g' ) } );
       else
          xticks( [ xmp xm x xp xpp ] );
          xticklabels( { num2str( xmp , '%1.1g' ) , 'x-' , 'x' , ...
                       'x+' , num2str( xpp , '%1.1g' ) } );
       end
       warning( 'on' , 'all' );
       hold off;
    end

    % check stopping condition- - - - - - - - - - - - - - - - - - - - - - -
    
    if abs( f1x ) <= eps
       break;
    end

    % restrict the interval based on sign of the derivative in xn - - - - -
   if f1x < 0
      xm = x;
      fxm = fx;
      f1xm = f1x;
   else
      xp = x;
      fxp = fx;
      f1xp = f1x;
   end

    if Interactive
       pause;
    end
    
    % iterate - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end

% end of main loop- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end  % the end- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




