%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Volume Integral of a Loader Quantity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function volint = volume_int(run_data, time_index, quantity)


	volint = run_data.rhomax(time_index)*trapz(linspace(0,1,length(run_data.qe(1,:))), quantity, 2);

end

