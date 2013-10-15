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

    %%%%%%%%%%%%%%%%%%%%%
    % Ions + Impurities %
    %%%%%%%%%%%%%%%%%%%%%

    Ai = param.compo.a
    Zi = param.compo.z
    
    %%%%%%%%
    % dens %
    %%%%%%%%

    ne = hello.prof.ne(t_index, :);
    gs2param.ne = interp1(x, ne, rho, 'spline');

    ni1 = hello.data.impur.impur(t_index, :, 1);
    gs2param.ni1 = interp1(x, ni1, rho, 'spline');
    
    ni2 = hello.data.impur.impur(t_index, :, 2);
    gs2param.ni2 = interp1(x, ni2, rho, 'spline');
    
    ni3 = hello.data.impur.impur(t_index, :, 3);
    gs2param.ni3 = interp1(x, ni3, rho, 'spline');
    
    ni4 = hello.data.impur.impur(t_index, :, 4);
    gs2param.ni4 = interp1(x, ni4, rho, 'spline');
    
    ni5 = hello.data.impur.impur(t_index, :, 5);
    gs2param.ni5 = interp1(x, ni5, rho, 'spline');
    
    %%%%%%%%
    % temp %
    %%%%%%%%

    te = hello.data.prof.te(t_index, :);
    gs2param.te = interp1(x, te, rho, 'spline');

    ti = hello.data.prof.ti(t_index, :);
    gs2param.ti1 = interp1(x, ti, rho, 'spline');
    
    gs2param.ti2 = gs2param.ti3 = gs2param.ti4 = gs2param.ti5 = gs2param.ti1;
    
    %%%%%%%%
    % zeff %
    %%%%%%%%

    zeff = hello.data.prof.zeff(t_index, :);
    gs2param.zeff = interp1(x, zeff, rho, 'spline');    
    
end
