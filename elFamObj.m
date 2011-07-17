function f = elFamObj(crLambda, p, elike, deriv)
% Computes the negated GEL objective function or its derivative for the given
% crLambda value. (0 for EL, -1 for exponential tilting).
% Used internally by other routines.
%
  switch crLambda
    case 0
      f = elComputePi(p,elike,deriv);
    case -1
      f = etiltComputePi(p,elike,deriv);
    otherwise
      f = gelComputePi(p, elike, deriv);
  end
  
  
function f = elComputePi(p, elike, deriv)
% Compute the objective function or its derivatives for empirical
% likelihood

  switch deriv
    case 0
      f = -sum(log(p));
    case 1
      f = -1 ./ p;
    case 2
      f = p .^ -2;
  end
end;

function f = etiltComputePi(p, elike, deriv)
% Compute objective function for exponential tilting

  switch deriv
    case 0
      f = -sum(p .* log(p));
    case 1
      f = log(p);
    case 2
      f = 1 ./ p;
  end
end;

function f = gelComputePi(p, elike, deriv)
% Works for general member of GEL family as long as there's no
% divide-by-zero problem (e.g. EL or ET)

  l = crLambda;

  switch deriv
    case 0
      f = sum(p .^ -l);
    case 1
      f = -l * p .^ (-1-l);
    case 2
      f = -l * (-1-l) * p .^ (-2-l);
  end

  % Note: omit the 1/(l * (1+l)) scaling factor (only the sign matters)
  if -1 < l && 0 > l
    f  = -f;
  end
  f = f / numel(p);

end;
  
end