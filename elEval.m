function res = elEval(elike, method, theta, p, lambda)
% Evaluates the empirical likelihood function at the given point
%
% Usage:
%    res = elEval(elike, theta, p, lambda);
%
% "elike" should be a structure obtained from elSetup(...)
%
% "method" is 'EL','ET', or a CR constant parameter
%
% "lambda" contains the Lagrange multipliers, but may be omitted
%          if you don't want the Hessian.
%
% "res" contains the following fields:
%    "res.f"      --  objective function
%    "res.fgrad"  --  gradient of objective
%    "res.H"      --  Hessian of the Lagrangian (empty if lambda missing)
%    "res.c"      --  constraints ( sum(pi(i) * m(i) )
%    "res.cJ"     --  Jacobian of moments per observation
%    "res.mom"    --  Moment conditions for each observation:
%                       an (nObs x nMom) matrix with one column for each
%                       moment condition and one row for each observation.
%

if (nargin <= 4)
  res = elSolve(elike, method, theta, p, 'eval');
else
  res = elSolve(elike, method, theta, p, 'eval', lambda);
end