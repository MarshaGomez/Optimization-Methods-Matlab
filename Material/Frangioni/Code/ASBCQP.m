function [ v , varargout ] =  ASBCQP( BCQP , varargin )

%function [ v , x , status ] = ASBCQP( BCQP , MaxIter )
%
% Apply the Active Set Method to the convex Box-Constrained Quadratic
% program
%
%  (P) min { (1/2) x^T * Q * x + q * x : 0 <= x <= u }
%
% encoded in the structure BCQP, where Q must be strictly positive
% definite.
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
 Author: Antonio Frangioni, Federico Poloni
 Date: 29-11-21
 Version 2.00
 Copyright Antonio Frangioni, Federico Poloni
 =======================================
%}

Interactive = false;  % if we pause at every iteration

Verbose = true;      % if we log every iteration

UpdateR = true;     % if we keep a Cholesky factorization of Q_{AA}
                     % updated rather than re-computing it every time

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

if UpdateR
   [ R , p ] = chol( BCQP.Q );
else
   [ ~ , p ] = chol( BCQP.Q );
end
if p > 0
   error( 'BCQP.Q not positive definite, this is not supported (yet)' );
end

if ~isempty( varargin )
   MaxIter = round( varargin{ 1 } );
   if ~ isscalar( MaxIter )
      error( 'MaxIter is not an integer scalar' );
   end
else
   MaxIter = 1000;
end

% initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

x = BCQP.u / 2;  % start from the middle of the box
v = 0.5 * x' * BCQP.Q * x + BCQP.q' * x;

% Because all constraints are box ones, the active set is logically
% partitioned onto the set of lower and upper bound constraints that are
% active, L and U respectively. Of course, L and U have to be disjoint.
% Since we start from the middle of the box, both the initial active sets
% are empty
L = false( n , 1 );
U = false( n , 1 );

% the set of "active variables", those that do *not* belong to any of the
% two active sets and therefore are "free", is therefore the complement to
% 1 : n of L union U; since L and U are empty now, A = 1 : n
A = true( n , 1 );

fprintf( 'Active Set method\n');
if Verbose
   fprintf( 'iter\tf(x)\t\t| B |\tI/O\n\n' );
end

i = 1;

% main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if ~ Interactive
   tic;
end

while true

    % output statistics - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if Verbose
       fprintf( '%4d\t%1.8e\t%d\t' , i , v , sum( L ) + sum( U ) );
    end
   
    % stopping criteria - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    if i > MaxIter
       status = 'stopped';
       break;
    end

    % solve the *unconstrained* problem restricted to A - - - - - - - - - -
    % the problem reads
    %
    %  min { (1/2) x_A' * Q_{AA} * x_A + ( q_A + u_U' * Q_{UA} ) * x_A }
    %    [ + (1/2) x_U' * Q_{UU} * x_U ]
    %
    % and therefore the optimal solution is
    %
    %   x_A^* = - Q_{AA}^{-1} ( q_A + u_U' * Q_{UA} )
    %
    % not that this actually is a *constrained* problem subject to equality
    % constraints, but in our case equality constraints just fix variables
    % (and anyway, any QP problem with equality constraints reduces to an
    % unconstrained one)
    
    xs = zeros( n , 1 );
    xs( U ) = BCQP.u( U );
    rhs = -( BCQP.q( A ) + BCQP.Q( A , U ) * BCQP.u( U ) );

    if UpdateR
       xs( A ) = R \ (R' \ rhs);
    else
       opts.SYM = true;     % tell Matlab Q_{AA} is symmetric positive
       opts.POSDEF = true;  % definite, it'll probably do the right thing
                            % and use Cholesky to solve the system
       xs( A ) = linsolve( BCQP.Q( A , A ) , rhs , opts );
    end

    if all( xs( A ) <= BCQP.u( A ) + 1e-12 & xs( A ) >= - 1e-12 )
       % the solution of the unconstrained problem is actually feasible

       % move the current point right there
       x = xs;

       % compute function value and gradient
       v = 0.5 * x' * BCQP.Q * x + BCQP.q' * x;
       g = BCQP.Q * x + BCQP.q;

       h = find( L & g < -1e-12 );
       if ~ isempty( h )
          uppr = false;          
       else
          h = find( U & g > 1e-12 );
          uppr = true;
       end
          
       if isempty( h )
          status = 'optimal';
          if Verbose
             fprintf( '\n' );            
          end
          break;
       else
          h = h( 1 , 1 );  % that's probably Bland's anti-cycle rule
          if UpdateR
              x1 = BCQP.Q(A, h);
              A( h ) = true;
              index_in_A = sum(A(1:h));
              R = cholinsert(R, x1, BCQP.Q(h,h), index_in_A);
              % norm(R'*R - BCQP.Q(A,A), 'fro') / norm(BCQP.Q(A,A), 'fro')  % the update ensures that this is small
          else
              A( h ) = true;        
          end

          if uppr
             U( h ) = false;
             if Verbose
                fprintf( 'O %d(U)\n' , h );            
             end
          else
             L( h ) = false;
             if Verbose
                fprintf( 'O %d(L)\n' , h );            
             end
          end
       end   
    else  % the solution of the unconstrained problem is *not* feasible
       % this means that d = xs - x is a descent direction, use it
       % of course, only the "free" part really needs to be computed

       d = zeros( n , 1 );
       d( A ) = xs( A ) - x( A );

       % first, compute the maximum feasible stepsize maxt such that
       %
       %   0 <= x( i ) + maxt * d( i ) <= u( i )   for all i

       ind = A & d > 0;  % positive gradient entries
       maxt = min( ( BCQP.u( ind ) - x( ind ) ) ./ d( ind ) );
       ind = A & d < 0;  % negative gradient entries
       maxt = min( [ maxt min( - x( ind ) ./ d( ind ) ) ] );
    
       % it is useless to compute the optimal t, because we know already
       % that it is 1, whereas maxt necessarily is < 1
       x = x + maxt * d;

       % compute function value
       v = 0.5 * x' * BCQP.Q * x + BCQP.q' * x;

       % update the active set(s)
       nL = A & x <= 1e-12;
       nU = A & x >= BCQP.u - 1e-12;
       
       if UpdateR
           indices = nL | nU;
           indices_in_A = indices(A); % restricts to columns of A
           R = choldelete(R, indices_in_A);
       end
       
       L( nL ) = true;
       A( nL ) = false;
       
       U( nU ) = true;
       A( nU ) = false;

       % norm(R'*R - BCQP.Q(A,A), 'fro') / norm(BCQP.Q(A,A), 'fro')  % the update ensures that this is small

       
       if Verbose
          fprintf( 'I %d+%d\n' , sum( nL ) , sum( nU ) );
       end
    end

    % iterate - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    i = i + 1;

    if Interactive
       pause;
    end
end

% end of main loop- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

fprintf( 'stop: %d iter, status = %s, fbest = %1.8e\n' , ...
         i , status , v );
if ~ Interactive
   toc
end

if nargout > 1
   varargout{ 1 } = x;
end

if nargout > 2
   varargout{ 2 } = status;
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function R = choldelete( R , k )
% update a Cholesky factorization A =R'*R after removing row+cols k from A
%
% k may be a logical vector or a sequence of indices

nr = length(R);
Q = eye(nr);

% convert k to a sequence of increasing indices
if islogical(k)
    k = find(k);
else
    k = sort(k);
end

for j = 1:length(k)
    % note that the columns in k change indices as we remove previous
    % columns
    [Q, R] = qrdelete(Q, R, k(j) - (j-1), 'col');
end
%resize
nr = size(R,2);
R = R(1:nr, :);

end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function R = cholinsert( R , x1 , x2 , k )
% update a Cholesky factorization A = R'*R after inserting a (single)
% row+col in position k
%
% Returns a Cholesky factorization of B, where B(k^c,k) = x1, B(k,k)=x2,
% B(k^c,k^c) = A

nr = length(R);

x1 = x1(:);

% first construct a quasi-triangular factorization
w = R' \ x1;
t = x2 - w'*w;
if t < 0
    error('Not numerically positive (semi)definite');
end
R = [R(:,1:k-1) w R(:,k:end); zeros(1,k-1) sqrt(t) zeros(1,nr-k+1)];

% then restore the triangular form via Givens transformations
for j = nr+1:-1:k+1
    [G, y] = planerot(R(j-1:j,k));
    R(j-1:j,k) = y;
    R([j-1,j],j:end) = G*R([j-1,j],j:end);
end

end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end  % the end- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




