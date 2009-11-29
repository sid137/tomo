function z = Wigner(n,rho,x,y)
% Wigner function
z=trace(Wmatrix(n,x,y)*(rho.'));

function wigm = Wmatrix(n,x,y)
% Wigner matrices
v=0;
for k=1:n+1
    for m=1:k
        v(k,m)=((-1)^(m-1))*sqrt((factorial(m-1))/(factorial(k-1)))*(sqrt(2.0)*(x-i*y))^(k-m)*LaguerreL(m-1,k-m,2*(x^2+y^2));
    end
end
u=tril(v)+tril(v,-1)';
wigm=(exp(-(x^2+y^2))/(pi))*u;

function w = LaguerreL(n, alpha, x)
%Fast Laguerre function, but limited to n=10.
if n==0
    w=1;
end;
if n==1
    w=alpha+1-x;
end;
if n==2
    w=(1/2)*(2+3*alpha+alpha^2-4*x-2*alpha*x+x^2);
end;
if n==3
    w=(1/6)*(6 + 11*alpha + 6*alpha^2 + alpha^3 - 18*x - 15*alpha*x - 3*x*alpha^2 + 9*x^2 + 3*alpha*x^2 - x^3);
end;
if n==4
    w=(1/24)*(24 + 50*alpha + 35*alpha^2 + 10*alpha^3 + alpha^4 - 96*x - 104*alpha*x - 36*x*alpha^2 - 4*x*alpha^3 + 72*x^2 + 42*alpha*x^2 + 6*(alpha^2)*x^2 - 16*x^3 - 4*alpha*x^3 + x^4);
end;
if n==5
    w=(1/120)*(120 + 274*alpha + 225*alpha^2 + 85*alpha^3 + 15*alpha^4 + alpha^5 - 600*x - 770*alpha*x - 355*(alpha^2)*x - 70*(alpha^3)*x - 5*(alpha^4)*x + 600*x^2 + 470*alpha*x^2 + 120*(alpha^2)*x^2 + 10*(alpha^3)*x^2 - 200*x^3 - 90*alpha*x^3 - 10*alpha^2*x^3 + 25*x^4 + 5*alpha*x^4 - x^5);
end;
if n==6
    w=(1/720)*(720 + 1764*alpha + 1624*alpha^2 + 735*alpha^3 + 175*alpha^4 + 21*alpha^5 + alpha^6 - 4320*x - 6264*alpha*x - 3480*(alpha^2)*x - 930*(alpha^3)*x - 120*(alpha^4)*x - 6*(alpha^5)*x + 5400*x^2 + 5130*alpha*x^2 + 1785*(alpha^2)*x^2 + 270*(alpha^3)*x^2 + 15*(alpha^4)*x^2 - 2400*x^3 - 1480*alpha*x^3 - 300*(alpha^2)*x^3 - 20*(alpha^3)*x^3 + 450*x^4 + 165*alpha*x^4 + 15*(alpha^2)*x^4 - 36*x^5 - 6*alpha*x^5 + x^6);
end; % valeurs verifiees et testees jusqu'a n=6 inclus
if n==7
    w=(1/5040)*(5040 + 13068*alpha + 13132*alpha^2 + 6769*alpha^3 + 1960*alpha^4 + 322*alpha^5 + 28*alpha^6 + alpha^7 - 35280*x - 56196*alpha*x - 35728*(alpha^2)*x - 11655*(alpha^3)*x - 2065*(alpha^4)*x - 189*(alpha^5)*x - 7*(alpha^6)*x + 52920*x^2 + 57834*(alpha)*x^2 + 24675*(alpha^2)*x^2 + 5145*(alpha^3)*x^2 + 525*(alpha^4)*x^2 + 21*(alpha^5)*x^2 - 29400*x^3 - 22330*alpha*x^3 - 6265*(alpha^2)*x^3 - 770*(alpha^3)*x^3 - 35*(alpha^4)*x^3 + 7350*x^4 + 3745*alpha*x^4 + 630*(alpha^2)*x^4 + 35*(alpha^3)*x^4 - 882*x^5 - 273*alpha*x^5 - 21*(alpha^2)*x^5 + 49*x^6 + 7*alpha*x^6 - x^7);
end;
if n==8
    w=(1/40320)*(40320 + 109584*alpha + 118124*alpha^2 + 67284*alpha^3 + 22449*alpha^4 + 4536*alpha^5 + 546*alpha^6 + 36*alpha^7 + alpha^8 - 322560*x - 554112*alpha*x - 390880*(alpha^2)*x - 147392*(alpha^3)*x - 32200*(alpha^4)*x - 4088*(alpha^5)*x - 280*(alpha^6)*x - 8*(alpha^7)*x + 564480*x^2 + 687456*alpha*x^2 + 340312*(alpha^2)*x^2 + 87780*(alpha^3)*x^2 + 12460*(alpha^4)*x^2 + 924*(alpha^5)*x^2 + 28*(alpha^6)*x^2 - 376320*x^3 - 332864*alpha*x^3 - 115920*(alpha^2)*x^3 - 19880*(alpha^3)*x^3 - 1680*(alpha^4)*x^3 - 56*(alpha^5)*x^3 + 117600*x^4 + 74620*alpha*x^4 + 17570*(alpha^2)*x^4 + 1820*(alpha^3)*x^4 + 70*(alpha^4)*x^4 - 18816*x^5 - 8176*alpha*x^5 - 1176*(alpha^2)*x^5 - 56*(alpha^3)*x^5 + 1568*x^6 + 420*alpha*x^6 + 28*(alpha^2)*x^6 - 64*x^7 - 8*alpha*x^7 + x^8);
end;
if n==9
    w=(1/362880)*(362880 + 1026576*alpha + 1172700*alpha^2 + 723680*alpha^3 + 269325*alpha^4 + 63273*alpha^5 + 9450*alpha^6 + 870*alpha^7 + 45*alpha^8 + alpha^9 - 3265920*x - 5973264*alpha*x - 4581036*(alpha^2)*x - 1932084*(alpha^3)*x - 491841*(alpha^4)*x - 77616*(alpha^5)*x - 7434*(alpha^6)*x - 396*(alpha^7)*x - 9*(alpha^8)*x + 6531840*x^2 + 8680608*alpha*x^2 + 4821768*(alpha^2)*x^2 + 1453284*(alpha^3)*x^2 + 257040*(alpha^4)*x^2 + 26712*(alpha^5)*x^2 + 1512*(alpha^6)*x^2 + 36*(alpha^7)*x^2 - 5080320*x^3 - 5058144*(alpha)*x^3 - 2064216*(alpha^2)*x^3 - 442260*(alpha^3)*x^3 - 52500*(alpha^4)*x^3 - 3276*(alpha^5)*x^3 - 84*(alpha^6)*x^3 + 1905120*x^4 + 1420524*alpha*x^4 + 418950*(alpha^2)*x^4 + 61110*(alpha^3)*x^4 + 4410*(alpha^4)*x^4 + 126*(alpha^5)*x^4 - 381024*x^5 - 207900*alpha*x^5 - 42210*(alpha^2)*x^5 - 3780*(alpha^3)*x^5 - 126*(alpha^4)*x^5 + 42336*x^6 + 16044*alpha*x^6 + 2016*(alpha^2)*x^6 + 84*(alpha^3)*x^6 - 2592*x^7 - 612*alpha*x^7 - 36*(alpha^2)*x^7 + 81*x^8 + 9*alpha*x^8 - x^9);
end;
if n==10
    w=(1/3628800)*(3628800 + 10628640*alpha + 12753576*alpha^2 + 8409500*alpha^3 + 3416930*alpha^4 + 902055*alpha^5 + 157773*alpha^6 + 18150*alpha^7 + 1320*alpha^8 + 55*alpha^9 + alpha^10 - 36288000*x - 69998400*alpha*x - 57537360*(alpha^2)*x - 26557640*(alpha^3)*x - 7611660*(alpha^4)*x - 1408890*(alpha^5)*x - 168840*(alpha^6)*x - 12660*(alpha^7)*x - 540*(alpha^8)*x - 10*(alpha^9)*x + 81648000*x^2 + 116672400*alpha*x^2 + 71122860*(alpha^2)*x^2 + 24193260*(alpha^3)*x^2 + 5029605*(alpha^4)*x^2 + 655200*(alpha^5)*x^2 + 52290*(alpha^6)*x^2 + 2340*(alpha^7)*x^2 + 45*(alpha^8)*x^2 - 72576000*x^3 - 79516800*alpha*x^3 - 36714720*(alpha^2)*x^3 - 9266880*(alpha^3)*x^3 - 1381800*(alpha^4)*x^3 - 121800*(alpha^5)*x^3 - 5880*(alpha^6)*x^3 - 120*(alpha^7)*x^3 + 31752000*x^4 + 26850600*alpha*x^4 + 9350040*(alpha^2)*x^4 + 1716750*(alpha^3)*x^4 + 175350*(alpha^4)*x^4 + 9450*(alpha^5)*x^4 + 210*(alpha^6)*x^4 - 7620480*x^5 - 4920048*alpha*x^5 - 1260000*(alpha^2)*x^5 - 160020*(alpha^3)*x^5 - 10080*(alpha^4)*x^5 - 252*(alpha^5)*x^5 + 1058400*x^6 + 506940*alpha*x^6 + 90510*(alpha^2)*x^6 + 7140*(alpha^3)*x^6 + 210*(alpha^4)*x^6 - 86400*x^7 - 29040*alpha*x^7 - 3240*(alpha^2)*x^7 - 120*(alpha^3)*x^7 + 4050*x^8 + 855*alpha*x^8 + 45*(alpha^2)*x^8 - 100*x^9 - 10*alpha*x^9 + x^10);
end;