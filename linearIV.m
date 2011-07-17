function res = linearIV(X, Z, y, meth, endogVars)
% Solves an instrumental variables problem using linear methods
%
% Can use either two-stage least squares (TSLS) or limited-information
% maximum likelihood (LIML). OLS is also an option, though it ignores Z
% may be inconsistent.
%
% Usage:
%
%   res = linearIV(X, Z, y, method, endogVars)
%
%   'X' is a matrix containing the structural variables
%
%   'Z' is a matrix containing the instruments. If some structural
%       variables are exogenous and to be used as instruments, the should
%       appear in Z as well as X.
%
%   'y' is a vector of the dependent variable
%
%   'method' is 'TSLS', 'LIML', or 'OLS'
%       (OLS is included for comparison, despite being inconsistent)
%
%   'endogVars' is a vector with length equal to the number of structural
%       variables. An entry contains 1 (or true) if the corresponding
%       column of X represents an endogenous variable.
%       endogVars is only used for LIML -- it may be omitted in other cases
%
% The returned structure 'res' contains:
%
%    'res.theta'   -- estimated coefficients on 'X'
%    'res.k'       -- in the LIML case, the value of the k parameter
%                     computed for the k-class estimator. Otherwise
%                     left blank.
%
  betaOLS = X \ y;
  
  % Two-stage least squares (consistent)
  if (strcmpi(meth,'tsls') || strcmpi(meth,'liml'))
    XpZ = X'*Z;
    PzX = XpZ * (Z \ X);
    PzXy = XpZ * (Z \ y);
    betaTSLS = PzX \ PzXy;
  end
  
  % LIML estimator via k-class approach
  if (strcmpi(meth,'liml'))
    if (~exist('endogVars') || length(endogVars) ~= size(X,2))
      error('Must provide endogVars vector with length = # struct vars');
    end
    endogVars = (endogVars == 1);   % convert to boolean
    
    % Compute the k to use for LIML.
    % Formulas based on Greene, 6th Ed, page 376
    
    % e1: residuals from projecting all endog onto ALL exog
    y0 = [y, X(:,endogVars)];
    e1 = y0 - Z*(Z \ y0);

    % e0: residuals from projecting all endog onto INCLUDED exog
    if (all(endogVars == true))
      % Special case with no included exogenous vars
      e0 = y0;
    else
      x0 = X(:, ~endogVars);
      e0 = y0 - x0*(x0 \ y0);
    end

    % Characteristic equation is the ratio of resid cross-prods
    w0 = e0'*e0;
    w1 = e1'*e1;
    detEq = w1 \ w0;

    % Get its smallest eigenvalue
    eigopts.disp = 0;
    k = eigs(detEq, 1, 'SM', eigopts);
    
    % Standard k-class formula
    liml2 = X'*y - k*X'*e1(:,1);
    liml1 = (1-k)*X'*X + k*PzX;
    betaLIML = liml1 \ liml2;

    res.k = k;
  end

  if (strcmpi(meth,'ols'))
    res.theta = betaOLS;
  elseif (strcmpi(meth,'tsls'))
    res.theta = betaTSLS;
  elseif (strcmpi(meth,'liml'))
    res.theta = betaLIML;
  else
    error('Unknown IV method');
  end
end