%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get approx index of an array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function appr_index = approx_index(array, value)

   appr_index = round(interp1(array,1:length(array),value));

end
