function s = char(p)
% CHAR - 
%   
if strcmp(p.typ,'sum')
  s = ['(' char(p.f1) ') + (' char(p.f2) ')'];
elseif strcmp(p.typ,'diff')
  s = ['(' char(p.f1) ') - (' char(p.f2) ')'];  
elseif strcmp(p.typ,'prod')
  s = ['(' char(p.f1) ') * (' char(p.f2) ')'];    
elseif strcmp(p.typ,'ratio')
  s = ['(' char(p.f1) ') / (' char(p.f2) ')'];    
elseif all(p.coeffs == 0)
  s = '0';
else
  if strcmp(p.typ,'polynomial')
    d = length(p.coeffs) - 1;
    s = [];
    nn = 0;
    for b = p.coeffs;
      cc(d+1-nn) = b;
      nn = nn + 1;
    end
    for a = cc;
      if a ~= 0;
	if ~isempty(s)
	  if a > 0
	    s = [s ' + '];
	  else
	    s = [s ' - '];
	    a = -a;
	  end
	end
	if a ~= 1 | d == 0
	  s = [s num2str(a)];
	  if d > 0
	    s = [s '*'];
	  end
	end
	if d >= 2
	  s = [s 'x^' int2str(d)];
	elseif d == 1
	  s = [s 'x'];
	end
      end
      d = d - 1;
    end
  elseif strcmp(p.typ,'gaussian')
    s = num2str(p.coeffs(1));
    s = ['Gaussian(' num2str(p.coeffs(1)) ',' ...
	 num2str(p.coeffs(2)) ',' ...
	 num2str(p.coeffs(3)) ')'];
  else
    s = ['*** char not yet implemented for' p.typ ' ***'];
  end    
end
