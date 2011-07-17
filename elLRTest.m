function [pval,lrstat] = elLRTest(elres, R, b, joint)
% Computes an empirical likelihood ratio test of the hypotheses "R*theta = b"
%
% Usage:
%
%    "testres = elLRTest(elres)"
% 
%    "testres = elLRTest(elres, R, b, joint)"
%
% Inputs:
%
%    "elres" is a result structure returned by elSolve
%
%    "R" is a K x nTheta matrix of linear hypotheses to test, while "b" is a
%        K x 1 vector of values to test R against. If "joint" is
%        omitted or false, then each row of the matrix is a separate
%        linear hypothesis to be tested of the form:
%
%            H0:  R(k,:) * theta = b
%
%        If "joint" is true, then the vector hypothesis: H0: hyp * theta = b
%        is tested.
%
%        If "R" is omitted, by default elLRTest tests separately that each
%        element of theta = 0. (So it is equivalent to "hyp = eye(nTheta)")
%
% Outputs:
%
%    "pval" is K x 1 vector for separate estimation or a scalar for joint
%        estimation. It is the p-value for the LR test: the probability under
%        the null that the LR statistic would be at least as large as the one
%        obtained. It comes from the CDF of a chi-squared(dof) distribution,
%        where "dof" is the rank of the hypothesis. If "joint" is false,
%        separate p-values are computed for each hypothesis.
%
%    "lrstat" the LR statistic obtained from each hypothesis:
%
%        LR = -2 * ( L(theta_restricted) - L(theta_unres) )
%
% Notes:
%
%    Each test requires re-solving the model, so it may take some time for a
%    large model.
%
  oldelike = elres.elike;
  l1 = elres.fval;
  if (numel(oldelike.linConsLHS) ~= 0)
    error(['elLRTest does not support models that already have ' ...
           'constraints']);
  end
  
  if nargin < 4
    % Default to separate tests
    joint = false;
  end
  if nargin < 2
    % Default hypotheses
    R = eye(oldelike.nTheta);
    b = zeros(size(R,1), 1);
  end
  if joint
    [pval, lrstat] = testone(R,b);
  else
    % Test each row separately
    nHyp = size(R,1);
    pval = zeros(nHyp, 1);
    lrstat = zeros(nHyp, 1);
    
    for iHyp=1:size(R,1)
      [p,h] = testone(R(iHyp,:), b(iHyp));
      pval(iHyp) = p;
      lrstat(iHyp) = h;
    end
  end
  
function [thisp,thislr] = testone(r1,b1)
  numCons = size(r1,1);
  elike = oldelike;
  elike.verbose = false;
  elike.linConsLHS = r1;
  elike.linConsRHS = b1;
  %  elike.solverOpts.initMu = 1e-6;
  elike.lambdaGuess = [elres.lagmult ; zeros(numCons,1)];
  thetaGuess = elres.theta;
  
  newres = elSolve(elike, 'EL', thetaGuess, elres.p);
  l2 = newres.fval;
  
  thislr = -2 * (l2 - l1);
  dof = rank(R);
  thisp = 1 - chi2cdf(thislr, dof);
end

end