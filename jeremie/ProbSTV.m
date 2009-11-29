function stv = ProbSTV(numPho,rS,nthS)
for n=0:numPho
stv(n+1,:,:)=(sqrt((nthS.^2).*(1+nthS).^2-(2.*nthS+1).^2.*(sinh(rS)).^2.*(cosh(rS)).^2)).^n ./((1+nthS).^2 + (2.*nthS+1).*(sinh(rS)).^2).^(n+(1/2)) .* LegendreL(n,(nthS.*(1+nthS))./(sqrt(nthS.^2*(1+nthS).^2-(2.*nthS+1).^2.*(sinh(rS)).^2.*(cosh(rS)).^2)));
end;

function leg = LegendreL(n,x)
if n==0
    leg=1;
end;
if n==1
    leg=x;
end;
if n==2
    leg=(1/2).*(-1 + 3.*x.^2);
end;
if n==3
    leg=(1/2).*(-3.*x + 5.*x.^3);
end;
if n==4
    leg=(1/8).*(3 - 30.*x.^2 + 35.*x.^4);
end;
if n==5
    leg=(1/8).*(15.*x - 70.*x.^3 + 63.*x.^5);
end;
if n==6
    leg=(1/16).*(-5 + 105.*x.^2 - 315.*x.^4 + 231.*x.^6);
end;
if n==7
    leg=(1/16).*(-35.*x + 315.*x.^3 - 693.*x.^5 + 429.*x.^7);
end;
if n==8
    leg=(1/128).*(35 - 1260.*x.^2 + 6930.*x.^4 - 12012.*x.^6 + 6435.*x.^8);
end;
if n==9
    leg=(1/128).*(315.*x - 4620.*x.^3 + 18018.*x.^5 - 25740.*x.^7 + 12155.*x.^9);
end;
if n==10
    leg=(1/256).*(-63 + 3465.*x.^2 - 30030.*x.^4 + 90090.*x.^6 - 109395.*x.^8 + 46189.*x.^10);
end;




%for n=0:numPho
%stv(n+1)=(sqrt(nthS^2*(1+nthS)^2-(2*nthS+1)^2*(sinh(rS))^2*(cosh(rS))^2))^n /((1+nthS)^2 + (2*nthS+1)*(sinh(rS))^2)^(n+(1/2)) * LegendreL(n,(nthS*(1+nthS))/(sqrt(nthS^2*(1+nthS)^2-(2*nthS+1)^2*(sinh(rS))^2*(cosh(rS))^2)));
%end;