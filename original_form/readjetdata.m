%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read JET PPFs directly using mdsplus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function jetppf = readjetdata(shotno, name_of_variable, public_or_private)
	switch name_of_variable
	case 'li'
		ddaname = 'EFIT';
		varname = 'XLI';
		ppfstring = 'PPF/MG2/XLI';
	case 'wdia'
		ppfstring = 'PPF/EFIT/WDIA'
	end  
;
	if nargin == 3
		private = public_or_private;
	else
		private = 0;
	end
	%[variable, time, radius];
	[variable, time, radius] = zmdsplusjet(shotno,ppfstring,private);
	jetppf.variable = variable;
	jetppf.time = time;
	jetppf.radius = radius;
end %
