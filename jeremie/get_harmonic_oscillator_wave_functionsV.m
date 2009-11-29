function wave_function=get_harmonic_oscillator_wave_functionsV(numPhotons, data)
x=data(:,2)';
coef = ((1/pi)^0.25)*exp(-(x.^2)/2);
for k=0:numPhotons
    wave_function(k+1,:)=coef.*exp(-i*k*data(:,1)').*psi(k,x)/sqrt((2^k)*factorial(k));
end


function w = psi(m,x)
% Defines wavefunctions up to order 10
if m==0
    w=1;
end;
if m==1
    w=2.*x;
end;
if m==2
    w=-2+4.*x.^2;
end;
if m==3
    w=-12.*x+8.*x.^3;
end;
if m==4
    w=12-48.*x.^2+16.*x.^4;
end;
if m==5
    w=120.*x-160.*x.^3+32.*x.^5;
end;
if m==6
    w=-120+720.*x.^2-480.*x.^4+64.*x.^6;
end;
if m==7
    w=-1680.*x+3360.*x.^3-1344.*x.^5+128.*x.^7;
end;
if m==8
    w=1680-13440.*x.^2+13440.*x.^4-3584.*x.^6+256.*x.^8;
end;
if m==9
    w=30240.*x-80640.*x.^3+48384.*x.^5-9216.*x.^7+512.*x.^9;
end;
if m==10
    w=-30240+302400.*x.^2-403200.*x.^4+161280.*x.^6-23040.*x.^8+1024.*x.^10;
end;