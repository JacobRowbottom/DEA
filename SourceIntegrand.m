function [f] = SourceIntegrand(s,xa,ya,xb,yb,Ox,Oy,Damp,ThetaLoc,ThetaStep)

% Source Integrand
% s from 0 to Dx
Dx=sqrt((xb-xa).^2+(yb-ya).^2);

xs=xa+s.*(xb-xa)./Dx;
ys=ya+s.*(yb-ya)./Dx;

R0=sqrt((Ox-xs).^2+(Oy-ys).^2);

SinTheta0=-((xb-xa).*(Ox-xs)+(yb-ya).*(Oy-ys))./(R0.*Dx);
CosTheta0=sqrt(1-SinTheta0.^2);
Theta0=asin(SinTheta0);

CosThetaCut=(abs(ThetaLoc-Theta0)<=0.5*ThetaStep).*CosTheta0;

f=exp(-Damp*R0).*CosThetaCut./(sqrt(Dx).*R0);

end

