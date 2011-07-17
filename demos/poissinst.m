function res = poissinst(X, y, Z, method, verbose)
% Solves a Poisson-IV model
%
% The endogeneity has the form of an omitted variable correlated with X:
%
%     E(y | X, eta) = exp(x*beta) * eta
%     E(eta | Z) = 1
%
% This is the model described in:
%
%  Mullahy, J. "Instrumental-Variable Estimation of Count Data Models:
%     Applications to Models of Cigarette Smoking." The Review of Economics and
%     Statistics (1997)
%
  [nObs,nX] = size(X);
  nZ = size(Z,2);
  if nargin < 5
    verbose = true;
  end
  
  % Obtain a starting guess from simple TSLS
  logY = log(y + 1) - exp(1);
  res1 = linearIV(X, Z, logY, 'tsls');
  theta0 = res1.theta;
  
  % Solve the model via the specified method, using Z instruments
  elike = elSetup(nObs, nZ, nX, @poissmom, 'verbose', verbose);
  resid = [];
  
  if strcmpi(method,'gmm')
    % Compute a first-stage weighting matrix based on the moments
    % evaluated at the initial guess if we're doing GMM
    nMom = nZ;
    mtmp = zeros(nObs,nMom);
    for im=1:nMom
      mtmp(:,im) = poissmom(theta0, elike, im, im==1);
    end
    elike.W1 = inv(mtmp'*mtmp / nObs);
  end
  
  res = elSolve(elike, method, theta0);

  % Note that the coefficient on the constant will not be exactly 1
  % due to the need to normalize E(eta) to 1.
  
function M = poissmom(theta, elike, jj, newTheta)
% Compute the Mullahy-form moment conditions
  if newTheta || true
    % Optimization: only recompute the residual when theta changes,
    % otherwise reuse the previously-computed value of 'resid'
    resid = y .* exp(-X * theta) - ones(nObs,1);
  end
  
  % Interact residual with instruments
  M = Z(:,jj) .* resid;
  
end;

end