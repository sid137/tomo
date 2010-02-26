function [density_matrix,likelihood_history]=maximum_likelihood_method2(phase, quadrature, photon_number, iteration)
    %
    % this program generates density matrix by means of maximum likelihood method
    %
    % phase: 1-d array of phase data
    % quadrature: 1-d array of quadrature data with same length as phase data
    % photon_number: number of photons for calculation of density matrix
    %
    % 2008/10/14 Y.Takeno
    %

tic;
    iteration = round(iteration);
    if iteration<10
    	iteration = 10;
    end

	likelihood_history = zeros(iteration, 1);
    density_matrix = single(normalize_density_matrix(eye(photon_number+1)));

    glob = get_global_setting();

    not_nan_index = find(not(isnan(phase) | phase==0));
    phase = phase(not_nan_index);
    quadrature = quadrature(not_nan_index);

    if size(phase, 1)>1
    	phase = phase';
    end
    if size(quadrature, 1)>1
    	quadrature = quadrature';
    end
    data = [phase(1,:); quadrature(1,:)];
    clear phase quadrature;

    data_length = size(data, 2);

    wave_function = get_harmonic_oscillator_wave_functions(photon_number, data);
    clear data;

    if data_length>300000
    	wave_function = single(wave_function);
    end

toc
    tic;
    for k=1:iteration

%slow
%{
        likelihood = 0;
        r_rho = zeros(photon_number+1, 'single');
        for d=1:data_length
            wf = wave_function(:,d);
            wfc = wf';
            pr = real(wfc*density_matrix*wf);
            r_rho = r_rho+wf*wfc/pr;
            likelihood = likelihood+log(pr);
        end
%}

%fast
%%{
        pr = real(dot(wave_function, density_matrix*wave_function));
        r_rho = (wave_function./kron(pr, ones(photon_number+1, 1)))*wave_function';
        likelihood = sum(log(pr));
%%}

        density_matrix = r_rho*density_matrix*(r_rho');
        density_matrix = density_matrix/trace(density_matrix);
        likelihood = likelihood/data_length;
        likelihood_history(k) = likelihood;
%        disp(['likelihood of ' num2str(k) 'th iteration: ' num2str(likelihood)]);
        if mod(k, 50)==0
            disp(['iteration: ' num2str(k)]);
            disp(['likelihood: ' num2str(likelihood, 16)]);
            toc;
        	pause(0.01);
        end
    end

    density_matrix = double(density_matrix);
end
