function opt = fieldopt(optstruct, optname, optdefault)
% Utility: get field value or default if field not present
  
  if isfield(optstruct, optname)
    opt = optstruct.(optname);
  else
    opt = optdefault;
  end
end