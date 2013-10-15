%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takes CRONOS profiles and calculates % 
%        GS2 input paramater           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initially will load full CRONOS results file. Will want to pass a structure
% All grids will be left on normalized toroidal flux. Miller will require reinterpolation.

% rho = normalized toroidal flux coordinate
% time = time in seconds
function gs2param = gs2_input(runnum, rho, time)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load full CRONOS results file %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    hello = load('~/run_results/JET/77933/run1_no_nbi/run1_no_nbi_resultat.mat');

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % x-grid and time index %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    x = hello.param.gene.x;
    t_index = approx_index(hello.data.gene.temps, time);

    %%%%%%%%
    % zeff %
    %%%%%%%%

    zeff = hello.data.prof.zeff(t_index, :);
    gs2param.zeff = interp1(x,zeff,rho)
    
end
