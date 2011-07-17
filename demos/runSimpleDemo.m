% runSimpleDemo.m --
%
% Generates data for and solves a simple instrumental variables model.
% Note that you should use elLinearIV if you have to solve such a model
% in practice (it's faster). This file is meant only as a simple example.
%

% Setup the location of the elike libraries
elLoad

% Generate the synthetic data from an linear model with endogenous X
randn('state', 1111); % keep the random seed stable
N = 1000;
L = 2;
J = 10;

% Parameters
truebeta = (1:L)';
truepi = ones(J,L-1);
sigeps = 0.5;
sigZ = 2;
sigX = 1;

% Exogenous and endogenous variables
% (delta introduces correlation between X and unobservable in y)
global X y Z;
Z = [ sigZ .* randn(N, J-1), ones(N, 1) ];
Xeps = sigX * randn(N,L-1);
X = [Z*truepi + Xeps, ones(N, 1)];
delta = 5.0 * ones(L-1,1);
y = X*truebeta + sigeps*randn(N,1) + Xeps*delta;

% Solve the model using standard automatic differentiation-based matElike
elike = elSetup(N, J, L, @demoMoments);
res = elSolve(elike, 'EL');

% Display the estimated result
names = strvcat('theta1', 'theta2');
elModelSumm(res, names);


% You can pass in an arbitrary CR lambda value instead of 'EL' or 'ET"
fprintf('Compare with CR lambda = -1.5\n');
elike = elSetup(N, J, L, @demoMoments);
res = elSolve(elike, -1.5);
elModelSumm(res, names);


% Test the null hypothesis that theta2 == 2
pval = elLRTest(res, [0,1], 2)

% Alternative method: specialized for linear instrumental vars
elivres = elLinearIV(X, Z, y, 'EL');
elivres.theta

