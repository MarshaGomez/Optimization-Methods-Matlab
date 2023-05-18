function v = optQf( Q , q )

%function v = optQf( Q , q )
%
% Given the data structure encoding a Quadratic function
%
%  f( x ) = (1/2) x^T * Q * x + q * x
%
% returns the value of f() in the "putative minimum" (or maximum), i.e.,
% the point
%
%   xStar = Q \ -q;
%
% Input:
%
% - Q ([ n x n ] real symmetric matrix, not necessarily positive
%   semidefinite): the quadratic part of f
%
% - q ([ n x 1 ] real column vector): the linear part of f
%
% Output:
%
% - f( xStar )
%
%{
 =======================================
 Author: Antonio Frangioni
 Date: 30-09-22
 Version 0.10
 Copyright Antonio Frangioni
 =======================================
%}

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

% "solve the problem" - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ignore if Q is singular

xStar = Q \ -q;

v = 0.5 * xStar' * Q * xStar + q' * xStar;

end