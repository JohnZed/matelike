function Mj = demoMoments(theta, elike, jj, newTheta)
% Simple IV model: interact Z with the computed epsilon for the given theta
  global X y Z;
  resid = y - X * theta;
  
  % Return a column vector with one element per observation
  Mj = Z(:,jj) .* resid;
end