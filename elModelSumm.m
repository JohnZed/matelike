function [stderrs,oidstat] = elModelSumm(res, names, quiet)
% Print some basic information about the EL model solution "res"
%
% Usage:
%
%   "[stderr,oidstat] = elModelSumm(res, names, [quiet])"
%
% Inputs:
%
%   "res" is a solution object returned by elSolve(...)
%   "names" is an optional vector of names for the elements of theta
%   "quiet" is an optional parameter -- if it's true, printing is suppressed
%
% Outputs:
%
%    "stderr" is a column vector of computed standard errors.
%    "oidstat" is the chi-squared statistic from an LR
%        test of overidentifying restrictions.
%
% Details:
%
%    Reports the estimated theta from the model, standard errors,
%    and a likelihood ratio-based test for the validity of the
%    overidentifying restrictions.
%
%    The standard errors are based on the normal asymptotic approximation
%    to the estimator of theta, using a variance matrix based on the moments
%    weighted by the computed "pi" values.
%
theta = res.theta;
elike = res.elike;
nTheta = elike.nTheta;

if ~exist('names') || isempty(names)
  names = num2str((1:elike.nTheta)');
end
if ~exist('quiet')
  quiet = false;
end
gmmMethod = false;
if strcmpi(res.meth, 'gmm')
  gmmMethod = true;
end

printNonquiet('nObs:           %d\n', elike.nObs);
printNonquiet('Moments:        %d\n', elike.nMom);
printNonquiet('Theta #elems:   %d\n', nTheta);
printNonquiet('Method:         %s\n', res.meth);
printNonquiet('\n');

% Test overidentifying restrictions via likelihood ratio test
oiddeg = elike.nMom - elike.nTheta;

if gmmMethod
  V = res.V;
  stderrs = sqrt(diag(V));
  if nargout > 1
    fprintf('Warning: OID test not implemented for GMM\n');
    oidstat = 0;
  end
else
  % Empirical likelihood LR-based over-id test
  % See Imbens (JBES 2002), Section 4.3
  p = res.p;
  Lunres = elike.nObs * log(1/elike.nObs);
  Lres = sum(log(p));
  oidLR = 2*(Lunres - Lres);
  chistat = chi2cdf(oidLR, oiddeg);

  if oiddeg
    printNonquiet('Degree of OverID:    %d\n', oiddeg);
    printNonquiet('Likelihood unres/restricted:    %12.5f / %12.5f\n',...
                  Lunres, Lres);
    printNonquiet('LR stat:         %9.5f\n', oidLR);
    printNonquiet('P-value:    %9.5f\n', 1 - chistat);
  else
    printNonquiet('Exactly-identified model\n');
    printNonquiet('Likelihood:      %9.5f\n', Lres);  
  end
  
  % Compute asymptotic approximation to standard errors
  fev = elEval(elike, res.meth, res.theta, res.p, res.lagmult);

  % This block of cJ is dM/dTheta (so it's already weighted by p)
  gamma = fev.cJ(1:elike.nMom,elike.nObs+1:end);

  % Weighted cross-product of moment conditions
  Bhat = (repmat(res.p,1,elike.nMom) .* fev.mom)' * fev.mom;

  % See Qin and Lawless (1994) pg 306
  V = inv(gamma' * (Bhat \ gamma)) / elike.nObs;
  stderrs = sqrt(diag(V));
  oidstat = oidLR;
end

printNonquiet('\nEstimated theta (and asymptotic stderr):\n\n');
for ii=1:elike.nTheta
  printNonquiet('%20s    %10.5f      (%8.5f)\n', ...
                names(ii,:), theta(ii), stderrs(ii));
end
printNonquiet('\n');

  
  function printNonquiet(varargin)
  % Prints only if 'quiet' is false
  if ~quiet
    fprintf(varargin{:});
  end
  end
end
