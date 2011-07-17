function optRes = dynamicPanel(nIndiv, nT, theta, method, usePsi2, varargin)
% Generate data for and solve a simple dynamic panel data model
% See "Imbens, 'Generalized Method of Moments and Empirical Likelihood',
% Journal of Business Economics and Statistics (2002)" for details.
%
% USAGE:
%
%       res = dynamicPanel(nIndiv, nT, theta, method, usePsi2)
%

MAXIT = 50;
PRINTLEV = 5;

% Use a burn-in period so that the results come from a long-term
% stationary distribution.
sigma_eta = 1;
sigma_eps = 0.3;
BURNIN = 200;
eta_i = sigma_eta * randn(nIndiv, 1);
Y = zeros(nIndiv, nT + BURNIN);
Y(:,1) = eta_i + sigma_eps * randn(nIndiv, 1);
for t=2:nT+BURNIN
  eps_t = sigma_eps * randn(nIndiv,1);
  Y(:,t) = theta*Y(:,t-1) + eta_i + eps_t;
end

% Get rid of burn-in periods
Y = Y(:,BURNIN+1:BURNIN+nT);

% Count the moment conditions
Nmom = (nT-1)*(nT-2) / 2;
if (usePsi2)
  Nmom = Nmom + (nT-2);
end
Ntheta = 1;

% Precompute the lags and diffs
Ylag = [repmat(NaN, nIndiv, 1), Y(:,1:end-1) ];
Ydiff = Y - Ylag;

% Store the time indices corresponding to each moment condition
% Each psi1 moment needs two time indices, stored in the columns of
% 'tidx', while each psi2 moment needs one time index, stored in tidx(:,1)
tidx = zeros(Nmom, 2);
ii = 1;
for t1=3:nT
  for t2=1:(t1-2)
    tidx(ii,1) = t1;
    tidx(ii,2) = t2;
    ii = ii + 1;
  end
end
if (usePsi2)
  for t1 = 3:nT
    tidx(ii,1) = t1;
    ii = ii+1;
  end
end


%
% Perform the actual estimation and time it
%
tic
  elike = elSetup(nIndiv, Nmom, Ntheta, @dpMoments, 'verbose', false);
  res = elSolve(elike, method, [],[]);
timeres = toc;

%
% Return the output (can be displayed with elModelSumm)
%
  optRes = res;

function M = dpMoments(theta, elike, mnum, newTheta)
%
% Evaluate moment condition 'mnum' for each observation
%

numPsi1 = (nT-1)*(nT-2) / 2;
if (mnum <= numPsi1)
  % Psi1 moment condition
  
  t1 = tidx(mnum,1);  t2 = tidx(mnum,2);
  
  lagdiff = Ydiff(:,t1) - theta*Ydiff(:,t1-1);
  M = Y(:,t2) .* lagdiff;
else
  % Psi2 moment condition

  t = tidx(mnum,1);
  M = Ydiff(:, t-1) .* ( Y(:,t) - theta*Y(:,t-1) );
end

end;


end