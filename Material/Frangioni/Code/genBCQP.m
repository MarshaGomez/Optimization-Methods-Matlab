function BCQP = genBCQP( n , varargin )

%function BCQP = genBCQP( n , actv , rank , ecc , seed , umin , umax )
%
% Produces a structure encoding a convex Box-Constrained Quadratic program
%
%  (P) min { (1/2) x^T * Q * x + q * x : 0 <= x <= u }
%
% Input:
%
% - n (integer, scalar): the size of the problem
%
% - actv (real, scalar, default 0.5): how many box constraints (as a
%   fraction of the number of variables n of the problems) the
%   unconstrained optimum will violate, and therefore we expect to be
%   active in the constrined optimum; note that there is no guarantee that
%   exactly acvt constraints will be active, they may be less or (more
%   likely) more, except when actv = 0 because then the unconstrained
%   optimum is surely feasible and therefore it will be the constrained
%   optimum as well
%
% - rank (real, scalar, default 1.1): Q will be obtained as Q = G^T G, with
%   G a m \times n random matrix with m = rank * n. If rank > 1 then Q can
%   be expected to be full-rank, if rank < 1 it will not
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
% - seed (integer, default 0): the seed for the random number generator
%
% - umin (real, scalar, default 8): the minimum value of each u_i
%
% - umax (real, scalar, default 12): the maximum value of each u_i
%
% Output: the BCQP structure, with the following fields:
%
% - BCQP.Q: n \times n symmetric positive semidefinite real matrix
%
% - BCQP.q: n \times 1 real vector
%
% - BCQP.u: n \times 1 real vector > 0
%
%{
 =======================================
 Author: Antonio Frangioni
 Date: 18-12-19
 Version 0.20
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
   actv = 0.5;
else
   actv = varargin{ 1 };
   if ~ isscalar( actv ) || ~ isreal( actv )
      error( 'actv not a real scalar' );
   end

   if actv < 0 || actv > 1
      error( 'actv must be in [0, 1]' );
   end
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
   ecc = varargin{ 3 };
   if ~ isscalar( ecc ) || ~ isreal( ecc )
      error( 'ecc not a real scalar' );
   end

   if ecc < 0 || ecc >= 1
      error( 'ecc must be in [0, 1)' );
   end
else
   ecc = 0.99;
end

if length( varargin ) > 3
   seed = varargin{ 4 };
   if ~ isscalar( seed )
      error( 'seed not a scalar' );
   end
   seed = round( seed );
else
   seed = 0;
end

rng( seed );

if length( varargin ) > 4
   umin = varargin{ 5 };
   if ~ isscalar( umin ) || ~ isreal( umin )
      error( 'umin not a real scalar' );
   end

   if umin <= 0
      error( 'umin must be > 0' );
   end
else
   umin = 8;
end

if length( varargin ) > 5
   umax = varargin{ 6 };
   if ~ isscalar( umax ) || ~ isreal( umax )
      error( 'umin not a real scalar' );
   end

   if umax <= umin
      error( 'umax must be > umin' );
   end
else
   umax = 12;
end

% generate u- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

BCQP.u = umin * ones( n , 1 ) + ( umax - umin ) * rand( n , 1 );

% generate Q- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

G = rand( round( rank * n ) , n );
Q = G' * G;

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

BCQP.Q = Q;

% generate q- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% We first generate the unconstrained minimum z of the problem in the form
%
%    min_x (1/2) ( x - z )^T * Q * ( x - z ) =
%          (1/2) x^T * Q * x - z^T * Q * x + (1/2) z^T * Q * z
%
% and then we set q = - z^T Q

z = zeros( n , 1 );

% outb( i ) = true if z( i ) will be out of the bounds
outb = rand( n , 1 ) <= actv;

% 50/50 chance of being left of lb or right of ub
lr = rand( n , 1 ) <= 0.5;
l = outb & lr;
r = outb & ~ lr;

% a random amount left of the lb (0)
z( l ) = - rand( sum( l ) , 1 ) .* BCQP.u( l );

% a random amount right of the ub (u)
z( r ) = BCQP.u( r ) .* ( 1 + rand( sum( r ) , 1 ) );

outb = ~ outb;  % entries that will be inside the bound
% pick at random in [ 0 , u ]
z( outb ) = rand( sum( outb ) , 1 ) .* BCQP.u( outb );

BCQP.q = - BCQP.Q * z;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end