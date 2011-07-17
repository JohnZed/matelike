function [profLL,profLR,profPV] = elValues(elres, thetaList, fixvars, dof, pIn, lambdaIn)
% Computes the profile empirical likelihood, p-value, or LR test statistic at
% various points. Useful for plotting profile likelihoods.
%
% The LR stat corresponds to the null hypothesis that:
%
%       "theta(fixvars) = thetaList(fixvars,i)"
%
% The test has "length(fixvars)" degrees of freedom unless otherwise specified.
%
% Usage:
%
%   "[ll,lr,pv] = elValues(elres, thetaList, fixvars, [dof])"
%
% Inputs:
%
%   "elres" is a result structure previously returned by elSolve(...)
%           solving the unconstrained problem.
%
%   "thetaList" is a matrix of size nFixed x nProfiles, where nFixed is the
%           length of "fixvars" and nProfiles is the number of points at
%           which to compute the profile likelihood.  Each column is a single
%           value of the parameter vector at which to compute the
%           likelihood.
%
%   "fixvars" is an integer vector of length nFixed indicating which elements
%           of theta should be held fixed. If fixvars is empty, it indicates
%           that all elements of theta are fixed (only the p's are
%           free). These all-theta-fixed problems can be solved faster.
%
%    "dof" specifies the number of degrees of freedom to use when computing
%        the p-values. (Only needed if p-values are computed)
%
% Outputs:
%
%    "ll" is a row vector containing the profile empirical log-likelihood at each
%          point in thetaList. (This is the log likelihood obtained after
%          maximizing over "p" and the theta elements not included in "fixvars",
%          but holding the other theta elements fixed.)
%
%    "lr" contains the profile empirical likelihood ratios (comparing the
%         restricted and unrestricted models)
%
%    "pv" contains the profile empirical p-value comparing the restricted and
%         unrestricted models(based on the chi-squared approximation to the
%         distribution of the LR stat)
%
  
  nProf = size(thetaList,2);
  profLL = zeros(1,nProf);
  elike = elres.elike;
  elike.verbose = false;

  nTheta = elike.nTheta;
  nObs = elike.nObs;
  nMom = elike.nMom;
  meth = elres.meth;
  prevTheta = elres.theta;
  prevLambda = elres.lagmult(:);
  prevP = elres.p;
  if exist('pIn'); prevP = pIn; end;
  if exist('lambdaIn'); prevLambda = lambdaIn; end;
  
  % Default to all thetas fixed
  if ~exist('fixvars') || isempty(fixvars)
    fixvars = 1:nTheta;
  else
    fixvars = int32(fixvars);
  end
  
  % Find the cressie-read lambda corresponding to method
  if isnumeric(meth)
    crLambda = meth;
  else
    switch upper(meth)
      case 'EL'
        crLambda = 0;
      case {'ET','ETILT'}
        crLambda = -1;
      case 'CUE'
        crLambda = -2;
      otherwise
        assert(false, 'Unknown method: %s\n', meth);
    end
  end
  
  if numel(fixvars) == nTheta
    % All theta variables are constrained, don't need to recompute
    % moment functions -- just do it once up front
    for ip=1:nProf
      theta = thetaList(:,ip);
      mom = zeros(nObs, nMom);
      for ii=1:nMom
        mom(:,ii) = elike.userCompMom(theta, elike, ii, (ii==1));
      end
      mom = mom';

      % Find the maximum empirical likelihood given this theta
      if crLambda ~= -2
        posvars = ones(nObs,1);
      else
        % CUE/Euclidean allow negative weights
        posvars = zeros(nObs,1);
      end
      opt = zipopts();
      opt.printstep = false;
      opt.maxWatchfail = 10;
        
      opt.maxIter = elike.maxIter;
      [x,info] = zipsolver(@reuseObj, @reuseCons, @reuseHess, prevP, opt, ...
                           prevLambda, posvars);
      if info.flag ~= 0
        fprintf('Computing profile likelihood element %d failed (flag=%d)\n', ...
                ip, info.flag);
      end
      
      % Store the result for this pass
      profLL(ip) = -reuseObj(x);
      prevP = x;
      prevLambda = info.lagmult(:);
    end
    
  else
    % Convert fixed variables to linear constraints
    % and do a full nonlinear optimization each time
    nFixed = length(fixvars);
    fixedmat = zeros(nFixed,nTheta);
    for ifix=1:nFixed
      fixedmat(ifix,fixvars(ifix)) = 1;
    end
    elike.linConsLHS = fixedmat;
    for ii=1:nProf
      thistheta = thetaList(:,ii);
      elike.linConsRHS = thistheta(fixvars);
      
      res = elSolve(elike, meth, prevTheta, prevP);
      profLL(ii) = res.fval;
      prevTheta = res.theta;
      prevP = res.p;
    end
  end

  % Possibly convert from profile log-likelihood to some other statistics
  ll_unc = elres.fval;

  if nargout > 1
    % compute likelihood ratio
    profLR = 2 * (ll_unc - profLL);
  end
  if nargout > 2
    % compute p-value from LR test
    profPV = 1 - chi2cdf(profLR, dof);
  end
  
function [f,g] = reuseObj(p)
% Computes the objective, re-using already-evaluated moments.
  f = elFamObj(crLambda, p, elike, 0);
  if nargout > 1
    g = elFamObj(crLambda, p, elike, 1);
  end
end

function [c,J] = reuseCons(p)
  c = [ mom * p; sum(p) - 1 ];
  J = [ mom; ones(1,nObs) ];
end

function H = reuseHess(p, objscale, lambda)
  eldiag = p .^ (-2);
  H = spdiags(objscale * eldiag, 0, nObs, nObs);
end
  
end