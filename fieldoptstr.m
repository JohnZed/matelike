function ismatch = fieldoptstr(optstruct, optname, tomatch)
  if ~isfield(optstruct, optname)
    ismatch = false;
  else
    ismatch = strcmpi(optstruct.(optname), tomatch);
  end
  
end