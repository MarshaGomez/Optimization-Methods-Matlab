function TF = oneDTestFunctions()

%function TF = oneDTestFunctions()
%
% Produces a cell array of function handlers, useful to test unconstrained
% one-dimensional optimization algorithms.
%
% Each function in the array has the following interface:
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
% - g (real, either scalar or a [ 1 x 2 ] matrix denoting an interval) is
%   the first optional argument. This also depends on x. if x == [] then
%   g is a [ 1 x 2 ] matrix denoting the standard interval in which f()
%   is supposed to be minimised (into which v is the minimum). If x ~= []
%   then g is a scalar containing the dervative g = f'(x) (or a
%   subgradient of f() in x if f() is not differentiable at x).
%
% - H (real, scalar) is the second optional argument. This must only be
%   specified if x ~= [], and it is the second derivative h = f''(x).
%   If no such information is available, the function throws error.
%
% The current list of functions is the following:
%
%  1 Taylored nasty 10-th degree polynomial with a local minimum in 0
%    (with value 0), a local minimum around -1, a saddle point around
%    +1, a local maximum around +2 and another local minimum around +3.
%    Thought to be plotted in [ -1.2 , 3.1 ], and minimised in
%    [ -0.5 , 3 ] (the global minimim there being 0). Note that, by
%    taking [ -1 , 3 ] instead, blind dycothomic search would hit jackpot
%    in two iterations by sheer luck.
%
%{
 =======================================
 Author: Antonio Frangioni
 Date: 08-11-18
 Version 1.01
 Copyright Antonio Frangioni
 =======================================
%}

TF = cell( 1 , 1 );
%TF{ 1 } = @(x) genericquad( [ 6 -2 ; -2 6 ] , [ 10 ; 5 ] , x );
TF{ 1 } = @custom;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% function [ v , varargout ] = genericquad( Q , q , x )
%  % generic quadratic function f(x) = x' * Q * x / 2 + q' * x
% 
% if isempty( x )  % informative call
%    if min( eig( Q ) ) > 1e-14
%       xStar = Q \ -q;
%       v = 0.5 * xStar' * Q * xStar + q' * xStar;
%    else
%       v = - Inf;
%    end
%    if nargout > 1
%       varargout{ 1 } = [ 0 ; 0 ];
%    end
% else
%    if ~ isequal( size( x ) , [ 2 1 ] )
%       error( 'genericquad: x is of wrong size' );
%    end
%    v = 0.5 * x' * Q * x + q' * x;  % f(x)
%    if nargout > 1
%       varargout{ 1 } = Q * x + q;  % \nabla f(x)
%       if nargout > 2
%          varargout{ 2 } = Q;       % \nabla^2 f(x)
%       end
%    end
% end
% end  % genericquad

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function [ v , varargout ] = custom( x )
% custom 10-th degree polynomial
% syms x
% f = @(x) 91 * x^2 / 30 - 19 * x^3 / 6 - 54 * x^4 / 25 ...
%   + 93 * x^5 / 23 - 23 * x^6 / 36 -  121 * x^7 / 93 ...
%   + 72 * x^8 / 91 - 13 * x^9 / 74 + 9 * x^10 / 640
%
% diff( f , x )
% 9*x^9)/64 - (117*x^8)/74 + (576*x^7)/91 - (847*x^6)/93
% - (23*x^5)/6 + (465*x^4)/23 - (216*x^3)/25 - (19*x^2)/2 + (91*x)/15
%
% diff( f , x , 2 )
% (81*x^8)/64 - (468*x^7)/37 + (576*x^6)/13 - (1694*x^5)/31
% - (115*x^4)/6 + (1860*x^3)/23 - (648*x^2)/25 - 19*x + 91/15

if isempty( x )  % informative call
   v = 0;
   if nargout > 1
      varargout{ 1 } = [ -0.5 3 ];
   end
else
   v = 91 * x^2 / 30 - 19 * x^3 / 6 - 54 * x^4 / 25 ...
     + 93 * x^5 / 23 - 23 * x^6 / 36 -  121 * x^7 / 93 ...
     + 72 * x^8 / 91 - 13 * x^9 / 74 + 9 * x^10 / 640;
   if nargout > 1
      g = 91 * x / 15 - 19 * x^2 / 2 - 216 * x^3 / 25 + 465 * x^4 / 23 ...
        - 23 * x^5 / 6 - 847 * x^6 / 93 + 576 * x^7 / 91 ...
        - 117 * x^8 / 74 + 9 * x^9 / 64;
      varargout{ 1 } = g;  % f'(x)
      if nargout > 2
         h = 91 / 15 - 19 * x - 648 * x^2 / 25 + 1860 * x^3 / 23 ...
           - 115 * x^4 / 6 - 1694 * x^5 / 31 + 576 * x^6 / 13 ...
           - 468 * x^7 / 37 + 81 * x^8 / 64;
         varargout{ 2 } = h;       % f''(x)
      end
   end
end
end  % custom

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end