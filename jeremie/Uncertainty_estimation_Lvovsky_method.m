%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program calculation uncertainty %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear variables and impose number format
clear
format long e;

numPhotons = 10; % max number of photons used for the calculation

xx=(-6:0.01:6); % creation axe des abscisses
datatemp(:,2)=xx'; % remplissage des data avec angle fixe et x variable
for m=451:500
for j=1:2000
    angletemp(1:1201)=angle(j+(m-1)*2000); % extraction d'une valeur fixe de angle, mesure j
    datatemp(:,1)=angletemp';
    wave_function = get_harmonic_oscillator_wave_functionsV(numPhotons,datatemp); % calcul des projecteurs Pi(theta_i, x) ou presque
    pr=real(dot(wave_function, rho*wave_function)); % calcul des probabilités p(theta_i,x), ça doit donner une gaussienne car proba gaussienne pour theta fixé et x qui varie !
    fresult=fit(xx',pr','gauss1'); % fit de p(theta_i,x) par une gaussienne
    deriv1=differentiate(fresult,xx'); % dérivée de la gaussienne par rapport à x
    sigma(j+(m-1)*2000)=sqrt(pr(701)./abs(deriv1(701))); % car en 701, ça correspond à x=1
end;
disp('Ecriture numéro  '), disp(m)
sigmatemp=sigma';
save('C:\sigmaLvovsky.txt','sigmatemp', '-ASCII', '-double');
end;


for numptstat=1:1
    disp('Iteration '),disp(numptstat)

    

% Clear variables and impose number format
clear wave_function R pr rho lnL valr sqzdB ropt probaPhotons nopt probaColPhotons minDiff antidB

numpts=1000000; % Nombre de points
posmin1arches=4; % Position du 1er min sur les arches de SQZ
posmin2arches=33; % Position du 2e min sur les arches de SQZ





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maximum likelihood program simplified for uncertainty %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numPhotons = 10; % max number of photons used for the calculation
% N.B: numPhotons is limited to 10 because I have entered manually the
% values of HERMITE(m,x) and Laguerre(m,alpha,x) for m = 0 to 10. It increases program speed by a
% factor of 4, useful when treating large data.
% Prompt the user for the desired number of bins
%disp('There are '), disp(numpt), disp(' points in the data.')
numBin=100;
% Associate an angle to each point
DeltaPhi=pi*numBin/(abs(posmin2arches-posmin1arches)*2*pi);
angle=(0:1:(numpts-1))-posmin1arches/100*numpts;
angle=angle'*DeltaPhi*2*pi/(numpts-1);
% Starting density matrix in the Fock basis
rho = eye(numPhotons+1)/(numPhotons+1);
% Number of iterations
numIt = 100;
% Computes the iterations
% lnL : likelihood logarithm
% rho is rotated at each iteration and lnL is diplayed
% when lnL is stationnary, the iterations can be stopped
data(:,1)=angle;
data(:,2)=normrnd(0,sigma_reconstruit);
% Calcul des numPhotons fonctions d'ondes une bonne fois pour toutes
wave_function = get_harmonic_oscillator_wave_functionsV(numPhotons,data);
clear data; % on n'a plus besoin des données car c'est dans wave_function maintenant
pr=real(dot(wave_function, rho*wave_function)); % calcul des probabilités p(theta_i,x_i) d'un seul coup !
lnL= sum(log(pr)); % calcul du log de vraisemblance
% disp(' ')
% disp('lnL initial = '), disp(lnL)
for m=1:numIt
    R=(wave_function./kron(pr, ones(numPhotons+1, 1)))*wave_function'; % calcul de la matrice R
    %sum(sum((rho-(R*rho*R/trace(R*rho*R)))/norm(rho)))
    rho=R*rho*R; % calcul de rho'
    rho=rho/trace(rho);
    pr=real(dot(wave_function, rho*wave_function)); % calcul des probabilités p(theta_i,x_i) avec le nouveau rho'
    %disp('lnL iteration number'), disp(m)
    %temp=lnL;
    lnL= sum(log(pr)); % calcul du log de vraisemblance avec le nouveau rho'
    %(lnL-temp)/temp
    %disp(lnL)
end;
clear pr m
% Some checks
% disp('Check that rho is hermitic, that is the value below (rho dagger - rho) is close to 0')
% disp('')
% normHrho=norm((rho'-rho),inf);
% disp(normHrho)
% pause
% disp('Check that rho''s eigenvalues (real part) are semi-definite positive')
% lambdarho=real(eig(rho));
% disp(lambdarho)
% pause
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
% plotx=(0:numPhotons);
% plot(plotx,probaPhotons,'*',plotx,ProbSTV(numPhotons,ropt,nopt), 'o')
% pause
statSQZ(numptstat,1)=sqzdB;
statSQZ(numptstat,2)=antidB;
save('C:\statSQZ.txt','statSQZ', '-ASCII', '-double');
end;