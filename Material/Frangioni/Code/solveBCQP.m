function [ v , varargout ] = solveBCQP( BCQP )

%function [ v , x ] = solveBCQP( BCQP )
%
% Solves the convex Box-Constrained Quadratic program
%
%  (P) min { (1/2) x^T * Q * x + q * x : 0 <= x <= u }
%
% encoded in the structure BCQP. The fields of the structure are:
%
% - BCQP.Q: n \times n symmetric positive semidefinite real matrix
%
% - BCQP.q: n \times 1 real vector
%
% - BCQP.u: n \times 1 real vector > 0
%
% Uses YALMIP and, in principle, whatever appropriate solver YALMIP has at
% its disposal for the task (parameters assume quadprog, though, so one
% could argue that using YALMPI is not tremendously useful).
%
% Input: the structure encoding the BCQP.
%
% Output:
%
% - v, the optimal value of the problem
%
% - optionally, x (n \times 1 real vector), the optimal solution
%
%{
 =======================================
 Author: Antonio Frangioni
 Date: 06-11-17
 Version 0.10
 Copyright Antonio Frangioni
 =======================================
%}

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

% construct YALMIP model- - - - - - - - - - - - - - - - - - - - - - - - - -

x = sdpvar( n , 1 );
F = [ 0 <= x <= BCQP.u ];
c = (1/2) * x' * BCQP.Q * x + x' * BCQP.q;

% solve BCQP- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% set parameters for cplex, your own solver may be different
%ops = sdpsettings( 'solver' , 'cplex' , 'cplex.threads' , 1 , ...
%                   'verbose' , 0 );
%ops = sdpsettings( 'solver' , 'cplex' , 'verbose' , 0 );
% set parameters for cbc
%ops = sdpsettings( 'solver' , 'cbc' , 'verbose' , 0 );
% last resort: MATLAB QP solver
ops = sdpsettings( 'solver' , 'QUADPROG' , 'verbose' , 0 );

optimize( F , c , ops );

v = value( c );

if nargout > 1
   varargout{ 1 } = value( x );
end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end