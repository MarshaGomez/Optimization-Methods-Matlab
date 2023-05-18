function [ Q , varargout ] = genQf( n , varargin )

%function [ Q , q , v ] = genQf( n , seed , rank , conv , ecc , dom , box )
%
% Produces the data structure encoding a Quadratic function
%
%  f( x ) = (1/2) x^T * Q * x + q * x
%
% Input:
%
% - n (integer, scalar): the size of the problem
%
% - seed (integer, default 0): the seed for the random number generator
%
% - rank (real, scalar, default 1.1): Q will be obtained as Q = G^T G, with
%   G a m \times n random matrix with m = rank * n (almost, see "conv"
%   below). If rank > 1 then Q can be expected to be full-rank, if
%   rank < 1 it will not
%
% - conv (real, scalar, default 1): a parameter controlling whether Q is
%   positive (semi)definite, negative (semi)definite, or indefinite. This
%   is done by actually computing Q as Q = G_+^T G_+ - G_-^T G_-. The m
%   rows of G are partitioned according to conv, i.e., m * conv rows are
%   put in G_+ while m * ( 1 - conv ) rows are put in G_-. Hence, conv = 1
%   means that Q is positive (semi)definite, conv = 1 means that Q is
%   negative (semi)definite, and any value in the middle means that Q is
%   indefinite, with "close to 1" meaning "almost positive (semi)definite",
%   and "close to 0" meaning "almost negative (semi)definite"
%
% - ecc (real, scalar, default 0.99): the eccentricity of Q, i.e., the
%   ratio ( \lambda_1 - \lambda_n ) / ( \lambda_1 + \lambda_n ), with
%   \lambda_1 the largest eigenvalue and \lambda_n the smallest one. Note
%   that this makes sense only if \lambda_n > 0, for otherwise the
%   eccentricity is always 1; hence, this setting is ignored if
%   \lambda_n = 0, i.e., Q is not full-rank (see above). An eccentricity of
%   0 means that all eigenvalues are equal, as eccentricity -> 1 the
%   largest eigenvalue gets larger and larger w.r.t. the smallest one
%
% - dom (real, scalar, default 1): possibly modifies the matrix by making
%   it more or less diagonally dominant. This is obtained by multiplying
%   each diagonal element by a random number uniformly picked in the
%   interval [ 2 * dom / 3 , 4 * dom / 3 ] (and therefore with average
%   value dom). If dom >> 1 this increases the diagonal dominance, if
%   dom << 1 this decreases it (if dom == 1, nothing is done)
%
% - box (real, scalar, default 1): this parameter controls the generation
%   of q. This is done by first randomly generating one (tentatively)
%   optimal solution x_* to the unconstrained minimization of f( x ) in the
%   symmetric hypercube of size 2 * abs( box ), i.e., by uniformly drawing
%   each of its entries in the interval [ - abs( box ) , abs( box ) ];
%   then, q is generated as - Q * x_*, so that Q * x_* + q = 0. This
%   ensures that f( x ) is bounded below if Q is positive semidefinite, and
%   bounded above if Q is negative semidefinite. However, we also want to
%   be able to create problems that are positive semidefinite but not
%   bounded below, which means that the system Q x + q = 0 must not have a
%   solution. This first of all requires Q to be low-rank (some of the
%   eigenvalues actually = 0), which can be obtained by putting rank < 1.
%   Then, if box < 0 the vector q obtained as above is modified by
%   multiplying each of its entries by a random number uniformly picked in
%   the interval [ 2 / 3 , 4 / 3 ], which makes it unlikely that the system
%   Q x = q still has a solution. However, this makes sense only if Q is
%   rank-deficient, hence this is done only if Q is not strictly positive
%   definite (for otherwise the system always has a solution).
%
% Output:
%
% - Q: n \times n symmetric real matrix, which is either positive
%      (semi)definite, negative (semi)definite, or indefinite according
%      to the value of the conv parameter
%
% - q: n \times 1 real vector, optional
%
% - v: real, optional: the optimal value of the optimizattion problem
%      min{ f( x ) : x \in \R^n }. This is -\infty if Q is not positive
%      semidefinite (see conv above). It is -\infty even if Q is positive
%      semidefinite but rank deficient (see rank above) and the system
%      Q x + q = 0 has no solution. If instead Q is positive semidefinite
%      (whatever its rank) and the system Q x + q = 0 has a solution x_*,
%      then x_* is an optimal solution to the problem and
%
%        v = f( x_* ) = (1/2) x_*^T * Q * x_* + q * x_*
%          = (1/2) [ x_*^T * Q * x_* + q * x_* ] + (1/2) q * x_*
%          = (1/2) x_*^T ( Q * x_* + q * ) + (1/2) q * x_*
%          = (1/2) q * x_*
%
%{
 =======================================
 Author: Antonio Frangioni
 Date: 15-08-22
 Version 0.10
 Copyright Antonio Frangioni
 =======================================
%}

if ~ isscalar( n ) || ~ isreal( n )
   error( 'n not a real scalar' );
end
n = round( n );
if n <= 0
   error( 'n must be > 0' );
end

if isempty( varargin )
   seed = 0;
else
   seed = varargin{ 1 };
   if ~ isscalar( seed ) || ~ isreal( seed )
      error( 'actv not a real scalar' );
   end

   seed = round( seed );
end

if length( varargin ) > 1
   rank = varargin{ 2 };
   if ~ isscalar( rank ) || ~ isreal( rank )
      error( 'rank not a real scalar' );
   end

   if rank <= 0
      error( 'rank must be > 0' );
   end
else
   rank = 1.1;
end

if length( varargin ) > 2
   conv = varargin{ 3 };
   if ~ isscalar( conv ) || ~ isreal( conv )
      error( 'conv not a real scalar' );
   end

   if ( conv < 0 ) || ( conv > 1 )
      error( 'conv must be in [0, 1)' );
   end
else
   conv = 1;
end

if length( varargin ) > 3
   ecc = varargin{ 4 };
   if ~ isscalar( ecc ) || ~ isreal( ecc )
      error( 'ecc not a real scalar' );
   end

   if ecc < 0 || ecc >= 1
      error( 'ecc must be in [0, 1)' );
   end
else
   ecc = 0.99;
end

if length( varargin ) > 4
   dom = varargin{ 5 };
   if ~ isscalar( dom )
      error( 'dom not a scalar' );
   end
   
   if dom < 0
      error( 'dom must be >= 0' );
   end
else
   dom = 1;
end

if length( varargin ) > 5
   box = varargin{ 6 };
   if ~ isscalar( box )
      error( 'box not a scalar' );
   end
   
   if box == 0
      error( 'box must not be 0' );
   end
else
   box = 1;
end

rng( seed );

% generate Q- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% first step: generate the appropriate rank and positive / negative
%             definiteness

r = round( rank * n );
p = round( r * conv );

if p > 0
   G = rand( p , n );
   Q = G' * G;
else
   Q = zeros( n , n );
end

if r > p
   G = rand( r - p , n );
   Q = Q - G' * G;
end

% second step: if dom ~= 1, modify the diagonal
% increase or decrease randomly each element by [ - 1/3 , 1/3 ] of its
% initial value, then multiply it by dom
if dom ~= 1
   D = diag( Q );
   D = D .* ( 1 + ( 2 * rand( n , 1 ) - 1 ) / 3 );
   D = dom * D;
   Q = spdiags( D , 0 , Q );
   Q = full( Q );
end

% compute eigenvalue decomposition
[ V , D ] = eig( Q );  % V * D * V' = Q , D( 1 , 1 ) = \lambda_n

if D( 1 , 1 ) > 1e-14
   % modify eccentricity only if \lambda_n > 0, for when \lambda_n = 0 the
   % eccentricity is 1 by default
   %
   % the formula is:
   %
   %                         \lambda_i - \lambda_n             2 * ecc
   % \lambda_i = \lambda_n + --------------------- * \lambda_n -------
   %                         \lambda_1 - \lambda_n             1 - ecc
   %
   % This leaves \lambda_n unchanged, and modifies all the other ones
   % proportionally so that
   %
   %   \lambda_1 - \lambda_n
   %   --------------------- = ecc      (exercise: check)
   %   \lambda_1 - \lambda_n

   d = diag( D );
   l = d( 1 ) * ones( n , 1 ) + ( d( 1 ) / ( d( n ) - d( 1 ) ) ) * ...
                                ( 2 * ecc / ( 1 - ecc ) ) * ...
                                ( d - d( 1 ) * ones( n , 1 ) );  

   Q = V * diag( l ) * V';
end

if nargout > 0
   % if so required generate q- - - - - - - - - - - - - - - - - - - - - - -
   %
   % we first generate the unconstrained minimum x_* of the problem in the
   % box [ - abs( box ) , abs( box ) ] and then we set q = - Q * x_*

   x = 2 * abs( box ) * rand( n , 1 ) - abs( box );

   q = - Q * x;

   % if so required, we now randomly destroy the alignment between q and
   % the image of Q so as to make it hard to solve Q x = q

   if ( box < 0 ) && ( D( 1 , 1 ) <= 1e-14 )
      q = q .* ( ( 4 / 3 ) * rand( n , 1 ) - 2 / 3 );
   end
   
   varargout{ 1 } = q;

   if nargout > 1
      % if so required compute v. - - - - - - - - - - - - - - - - - - - - -
      % v is finite-valued only if either Q is strictly positive definite
      % or it is positive semidefinite but q has been constructed in such
      % a was that Q * x + q = 0 has a solution (that is the x we have)

      if ( D( 1 , 1 ) > 1e-14 ) || ...
         ( ( D( 1 , 1 ) > -1e-14 ) && ( box > 0 ) )
         v = q' * x / 2;
      else
         v = - Inf;
      end
       
      varargout{ 2 } = v;
   end
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end