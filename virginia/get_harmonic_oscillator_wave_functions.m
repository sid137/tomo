function wave_function=get_harmonic_oscillator_wave_functions(photon_number, data)
    %
    % photon_number: number of photons
    % data: 2-d array (data(1,:): phase array, data(2,:): quadrature array)
    %
    % original script by H.Yonezawa-san
    % modified by Y.Takeno 2008/10/14
    %

    glob = get_global_setting();

    index = find(data(2,:)>8.0);
    data(2,index) = 8.0;
    index = find(data(2,:)<-8.0);
    data(2,index) = -8.0;

    x = sqrt(1/glob.h_bar)*data(2,:);
    hermite = ones(photon_number+1, length(x));
    hermite(2,:) = 2*x;
    coef = ((2/pi/glob.h_bar)^0.25)*exp(-x.*x/(2*glob.h_bar));
    wave_function = zeros(size(hermite), 'single');
    for pn=0:photon_number
        if pn>1
            hermite(pn+1,:) = 2*x.*hermite(pn,:)-2*(pn-1)*hermite(pn-1,:);
        end
        wave_function(pn+1,:) = coef.*exp(i*pn*data(1,:)).*hermite(pn+1,:)/sqrt(2^pn*factorial(pn));
    end
end
