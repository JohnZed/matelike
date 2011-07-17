function res = elLinearIV(X, Z, y, elMethod, varargin)
% Solves a linear instrumental variables problem with EL
%
% Usage:
%
%   "res = elLinearIV(X, Z, y, elMethod)"
%
% Inputs:
%
%   "X" is a matrix containing the structural variables
%   "Z" is a matrix containing the instruments
%   "y" is a vector of the dependent variable
%   "elMethod" specifies the method to use ('EL','ET', or a Cressie-Read lambda
%       value). You can omit this argument to use the default 'EL'
%
%  The returned structure "res" is the same as that returned by elSolve. The
%  field "res.theta" contains the estimates for the coefficients on X. Note that
%  it contains a field "elike" with the problem setup.  This may be useful to
%  pass to "elEval(...)", etc.
%

if (nargin < 4)
% Default to empirical likelihood
    elMethod = 'EL';
end


nMom = size(Z,2);
nObs = size(X, 1);
assert(nObs == size(Z,1), ...
       'Instrument and structural vars have different nObs');
nTheta = size(X,2);
resid = [];

% Set up the empirical likelihood problem
elike = elSetup(nObs, nMom, nTheta, @livMoments, varargin{:});

% Use custom implementations of jacobian and hessian, not automatic diff
elike.userCompJac = @livGrad;
elike.userCompHessTheta = @livHess;

spZ = spones(Z);
if nnz(spZ) > 0.25 * numel(Z)
  % Use dense matrices internally, since Z is fairly dense
  elike.momSparsePat = [];
else
  elike.momSparsePat = spones(Z);
end

% Use OLS estimate as staring guess
thetaOls = X \ y;

% Actually solve the EL problem
res = elSolve(elike, elMethod, thetaOls, [], varargin{:});

function M = livMoments(theta, elike, mnum, newTheta)
% Compute the moments, given theta

if (newTheta)
  resid = (y - X*theta);
end

M = resid .* Z(:,mnum);
end;

function grads = livGrad(theta, elike, mnum, newTheta)
% Compute the gradient of moment number 'mnum' for each observation
    
grads = -X .* repmat(Z(:,mnum),1,nTheta);
end;

function H = livHess(theta, elike, mnum, newTheta)
% Compute the (trivially 0) Hessian for moment condition 'mnum'

H = sparse(nTheta,nTheta);
end;


end