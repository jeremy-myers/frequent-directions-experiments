clear all;

%% Experiment parameters
k = 5;
l = 20;

%% A is the data matrix
A = readmatrix( 'matrix.txt', 'Delimiter', ' ' );
[Au,As,Av] = svd(A,'econ');

%% Error right hand sides (depend on data)
fnormAk     = norm( A - Au(:,1:k)*(As(1:k,1:k)*Av(:,1:k)'), 'Fro' );
coverr_rhs  = (1/(l-k));
projerr_rhs = (l/(l-k));

%% Load results
T = readmatrix( 'shrinkt.txt', 'Delimiter', ' ' );
F = readmatrix( 'shrinkf.txt', 'Delimiter', ' ' );

%% Need SVDs
[Tu,Ts,Tv] = svd( T, 'econ' );
[Fu,Fs,Fv] = svd( F, 'econ' );

%% Cov error
coverrT_lhs = norm( A'*A - T'*T );
coverrF_lhs = norm( A'*A - F'*F );

%% Proj error
projerrT_lhs = norm( A - A * (Tv*Tv') , 'Fro' )^2;
projerrF_lhs = norm( A - A * (Fv*Fv') , 'Fro' )^2;

%% Print results
% Covariance error
fprintf( '\nCov-err\n' );
fprintf( '  with shrinking:   ' );
fprintf( ' %f <= %f\n', coverrT_lhs/norm( A, 'Fro')^2, coverr_rhs );
fprintf( '  without shrinking:' );
fprintf( ' %f <= %f\n', coverrF_lhs/norm( A, 'Fro')^2, coverr_rhs );

% Projection error
fprintf( '\nProj-err\n' );
fprintf( '  with shrinking:   ' );
fprintf( ' %f <= %f\n', projerrT_lhs/fnormAk^2, projerr_rhs );
fprintf( '  without shrinking:' );
fprintf( ' %f <= %f\n', projerrF_lhs/fnormAk^2, projerr_rhs );

% Other residuals
Rt = zeros(k,1);
Rf = zeros(k,1);

for i = 1:k
    av = A*Tv(:,i);
    u  = 1/Ts(i,i) * av;
    Rt(i) = norm( [ av - Ts(i,i)*u; A'*u - Ts(i,i)*Tv(:,i) ] );
    
    av = A*Fv(:,i);
    u  = 1/Fs(i,i) * av;
    Rf(i) = norm( [ av - Fs(i,i)*u; A'*u - Fs(i,i)*Fv(:,i) ] );
end

fprintf('\nResidual norms\n' );
fprintf('  with shrinking:\n' );
fprintf('    %f\n', Rt/As(1,1) );
fprintf('\n')
fprintf('  without shrinking:\n' );
fprintf('    %f\n', Rf/As(1,1) );