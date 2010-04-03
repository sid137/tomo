% Clear variables and impose number format
clear
format long e;
% Loads special polynomials package
% maple('with','orthopoly');
disp('')
disp('*******************************')
disp('* Frequency filtering program *')
disp('*******************************')
disp('')
% Prompts user for the sampling frequency in MHz (Ms/s).
%Fs = 1000000*input('Sampling frequency (Ms/s) ?');
Fs = 1000000*100;

% Loads the squeezing data.
load data.txt

phase = data(:,2)*2*pi/360;

sqz = data(:,1);
number_of_points = size(sqz,1);
NfftPt=number_of_points; % Number of points used for fft
% Loads the shot noise data
shot = data(:,3); % N.B: data MUST be a number_of_points*1 array, that is data must be in one column !!!

% Starts the filtering program for the sqz
%time = 0:(1/Fs):((NfftPt-1)/Fs);
% plot(time,sqz)
% title('Signal photocurrent')
% xlabel('Time (s)')
% ylabel('Quadrature amplitude (a.u.)')
% Computes the FFT
FFTsqz = fft(sqz,NfftPt);
Psqz = FFTsqz.*conj(FFTsqz)/NfftPt;
%freq = (0:Fs/NfftPt:Fs/2);
disp('Data ranging from 0 to')
disp(Fs/(2000000))
disp('MHz')
disp(' ')
% Asks user for low cut and high cut frequencies
%lowCutFreq=1000*input('Low cut filter frequency ? (kHz)  ');
lowCutFreq=1000*100;
%highCutFreq=1000000*input('High cut filter frequency ? (MHz)  ');
highCutFreq=1000000*3;
lowCutPt=fix(lowCutFreq/(Fs/NfftPt));
highCutPt=fix(highCutFreq/(Fs/NfftPt));
% Replaces filtered frequencies by 0;
if (lowCutPt>0) && (lowCutPt<((NfftPt/2)+1))
    FFTsqz(1:lowCutPt)=0;
    FFTsqz(NfftPt-lowCutPt+1:NfftPt)=0;
end;
if (highCutPt>0) && (highCutPt<((NfftPt/2)+1))
    FFTsqz(highCutPt:(NfftPt-highCutPt))=0;
end;
% N.B: another way instead of putting 0 for the FFT values would be to
% replace them by the electronic noise from the card at those frequencies.
% That would be equal to putting filters just before the card.
% Calculates the filtered data
sqz=real(ifft(FFTsqz,NfftPt));
% Just a check, can be removed to go faster.
%FFTsqz=fft(sqz,NfftPt);
%Psqz = FFTsqz.*conj(FFTsqz)/NfftPt;
%plot(freq(lowCutPt+1:highCutPt),Psqz(lowCutPt+1:highCutPt))
%title('Frequency content of filtered signal')
%xlabel('Frequency (Hz)')
%ylabel('Power spectrum (a.u.)')
%clear FFTsqz time Psqz freq
%pause

% Starts the filtering program for the shot using the same frequencies for
% the filters
time = 0:(1/Fs):((NfftPt-1)/Fs);
% plot(time,shot)
% title('Signal photocurrent')
% xlabel('Time (s)')
% ylabel('Quadrature amplitude (a.u.)')
% Computes the FFT
FFTshot = fft(shot,NfftPt);
Pshot = FFTshot.*conj(FFTshot)/NfftPt;
freq = (0:Fs/NfftPt:Fs/2);
lowCutPt=fix(lowCutFreq/(Fs/NfftPt));
highCutPt=fix(highCutFreq/(Fs/NfftPt));
% Replaces filtered frequencies by 0;
if (lowCutPt>0) && (lowCutPt<((NfftPt/2)+1))
    FFTshot(1:lowCutPt)=0;
    FFTshot(NfftPt-lowCutPt+1:NfftPt)=0;
end;
if (highCutPt>0) && (highCutPt<((NfftPt/2)+1))
    FFTshot(highCutPt:(NfftPt-highCutPt))=0;
end;
% N.B: another way instead of putting 0 for the FFT values would be to
% replace them by the electronic noise from the card at those frequencies.
% That would be equal to putting filters just before the card.
% Calculates the filtered data
shot=real(ifft(FFTshot,NfftPt));
% Just a check, can be removed to go faster.
%FFTshot=fft(shot,NfftPt);
%Pshot = FFTshot.*conj(FFTshot)/NfftPt;
%plot(freq(lowCutPt+1:highCutPt),Pshot(lowCutPt+1:highCutPt))
%title('Frequency content of filtered shot')
%xlabel('Frequency (Hz)')
%ylabel('Power spectrum (a.u.)')
%clear FFTshot time Pshot freq lowCutPt highCutPt
%pause


disp('')
disp('******************************')
disp('* Maximum likelihood program *')
disp('******************************')
disp('')
numPhotons = 10; % max number of photons used for the calculation
% N.B: numPhotons is limited to 10 because I have entered manually the
% values of HERMITE(m,x) and Laguerre(m,alpha,x) for m = 0 to 10. It increases program speed by a
% factor of 4, useful when treating large data.
% Prompt the user for the desired number of bins
%disp('There are '), disp(number_of_points), disp(' points in the data.')
%numBin = input('Number of bins ? (multiple of 4) ');
numBin=100;
% Plots the shot noise
plot(shot)
%pause
% Plots the sqz
plot(sqz)
%pause
% The mean of the shot noise (corresponding to
% the mean of electronics noise) is substracted to the sqz. The SQL
% given by the standard deviation of the shot noise is used to normalize
% the data.
zero=mean(shot(:,1));
SQL=std(shot(:,1)-zero);
shotVar=var((shot-zero)/(SQL*sqrt(2.0)));
dataSQZ=((sqz-zero)/(SQL*sqrt(2.0)));
%clear sqz shot SQL
% Defines the marginal distribution. In order to do this, the vector of
% data is cut into numBin bins of approximately number_of_points/numBin points. The
% output is an array of size numBin x sizeBin in which each column is a
% marginal distribution. 
sizeBin = fix(number_of_points/numBin);
for j=1:numBin
   Pm(:,j)=(dataSQZ(sizeBin*(j-1)+1:sizeBin*j,1));
end
clear j
% For each of the marginal distribution (say for each of the columns 
% of the arrays of marginal distributions) it computes the variance. The 
% plot of the variances as a function of the bins is evaluated.
noise = zeros(1,numBin);
for k=1:numBin
    noise(k)=var(Pm(:,k));  
end
clear k
%plot(noise,'.')
%pause
plot(10*log10(noise/0.5))

% Starting density matrix in the Fock basis
rho = eye(numPhotons+1)/(numPhotons+1);
% Number of iterations
numIt = 100;
% Computes the iterations
% lnL : likelihood logarithm
% rho is rotated at each iteration and lnL is diplayed
% when lnL is stationnary, the iterations can be stopped

data(:,1)=phase;
data(:,2)=dataSQZ;
% Calcul des numPhotons fonctions d'ondes une bonne fois pour toutes
wave_function = get_harmonic_oscillator_wave_functionsV(numPhotons,data);
clear data; % on n'a plus besoin des donn�es car c'est dans wave_function maintenant
pr=real(dot(wave_function, rho*wave_function)); % calcul des probabilit�s p(theta_i,x_i) d'un seul coup !
lnL= sum(log(pr)); % calcul du log de vraisemblance
disp(' ')
disp('lnL initial = '), disp(lnL)
for m=1:numIt
    R=(wave_function./kron(pr, ones(numPhotons+1, 1)))*wave_function'; % calcul de la matrice R
    %sum(sum((rho-(R*rho*R/trace(R*rho*R)))/norm(rho)))
    rho=R*rho*R; % calcul de rho'
    rho=rho/trace(rho);
    pr=real(dot(wave_function, rho*wave_function)); % calcul des probabilit�s p(theta_i,x_i) avec le nouveau rho'
    disp('lnL iteration number'), disp(m)
    %temp=lnL;
    lnL= sum(log(pr)); % calcul du log de vraisemblance avec le nouveau rho'
    %(lnL-temp)/temp
    disp(lnL)
    if (mod(m,5) == 0)
        save lnL.mat lnL;
    end
end;
clear pr
% Some checks
disp('Check that rho is hermitic, that is the value below (rho dagger - rho) is close to 0')
disp('')
normHrho=norm((rho'-rho),inf);
disp(normHrho)
%pause
disp('Check that rho''s eigenvalues (real part) are semi-definite positive')
lambdarho=real(eig(rho));
disp(lambdarho)
%pause
% Results
% Computes the diagonal values of rho = probabilities of having 0, 1, 2, ...
% photons
probaPhotons=real(diag(rho));
% Computes the squeezing parameter r that best suits data, so that the 
% theoretical and experimental probability of having 2 photons agree.
% N.B: r goes from 0 to 0.6, that is from 0 dB of squeezing to -5.2 dB.
% For more squeezing, the calculation of r should be done differently or increase value of r.
% nth goes from 0 to 0.4 can be increased as well.
clear valr
valr(1,1:600)=(0:599)/1000;
valr(1,1)=0.00001;
valr=valr'; % column vectro with all values of r
probaColPhotons = kron(ones(1,600),probaPhotons); % set of 600 columns, each column having probaPhotons
v=0;
nopt=0;
ropt=0;
minDiff=1000000;
for nth=0:0.001:0.4
    v=ProbSTV(numPhotons,valr,nth);
    [minCol,indexCol]=min(sum((v-probaColPhotons).^2));
    if minCol < minDiff
        minDiff = minCol;
        ropt = valr(indexCol);
        nopt = nth;
    end;
end;
clear nth v minCol indexCol
% Computes the resulting squeezing and antisqz levels
disp('Squeezing found, in dB (negative if squeezing): ')
sqzdB= 10*log10((2*nopt+1)*exp(-2*ropt));
disp(sqzdB)
disp('Anti-squeezing found, in dB: ')
antidB= 10*log10((2*nopt+1)*exp(2*ropt));
disp(antidB)
% Plots the probability of having 0, 1, 2, ... photons, theory (o) and
% experiment (*)
plotx=(0:numPhotons);
plot(plotx,probaPhotons,'*',plotx,ProbSTV(numPhotons,ropt,nopt), 'o')
%pause
% 3-D Graph of the Wigner function
numPtGraph = input('Number of disivions on each axis of the Wigner graph ?   ');
maxXWig = input('Horizontal axes ranging from -a to a; a = ?   ');
[XWig,YWig] = meshgrid([-maxXWig:(2*maxXWig/numPtGraph):maxXWig]);
PlotWig = zeros(numPtGraph+1, numPtGraph+1);
for k1=1:(numPtGraph+1)
    for k2=1:(numPtGraph+1)
    PlotWig(k1,k2)=Wigner(numPhotons,rho,XWig(k1,k2),YWig(k1,k2));
    end
end
clear k1 k2
PlotWig=real(PlotWig);
disp('Graph results')
figure(1)
plot(dataSQZ)
%print('-f1', '-r600', '-depsc', 'E:\temp\Photocurrent');
figure(2)
surf(real(PlotWig))
shading interp
%print('-f2', '-r600', '-depsc', 'E:\temp\Wig3D');
figure(3)
contourf(PlotWig,50)
%print('-f3', '-r600', '-depsc', 'E:\temp\ProjWig');
absRho=abs(rho);
figure(4)
bar3(absRho,'detached')
%print('-f4', '-r600', '-depsc', 'E:\temp\HistRho');
figure(5)
contourf(XWig,YWig,PlotWig,[max(max(PlotWig))*0.882,max(max(PlotWig))*0.882])
axis([-1 1 -1 1])
set(gca,'xtick',XWig(1,:))
set(gca,'ytick',YWig(:,1))