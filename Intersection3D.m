function [Xstar, Ystar, Zstar, Angle] = Intersection3D(xai,yai,zai,xbi,ybi,zbi,xaj,yaj,zaj,xbj,ybj,zbj)

%Inputs: vertex coordinates (xai,yai) and (xbi,ybi) of edge i and (xaj,yaj) and (xbj,ybj) of edge j 
%Transition is j to i

%%Outputs:%%

% Angle(i,j)= internal angle in polygon between edges 1 and 2 
% (Xstar,Ystar) are coordinates of position where two edges meet at angle Angle
% (Xstar,Ystar) set by default to (0,0) if parallel edges 
% L0 = min length of perpendicular line between two edges

 Dxi=(xbi-xai);              % x-dist btwn vertices on edge i
 Dyi=(ybi-yai);              % y-dist btwn vertices on edge i
 Dzi=(zbi-zai);              % z-dist btwn vertices on edge i
 Dxj=(xbj-xaj);               % x-dist btwn vertices on edge j
 Dyj=(ybj-yaj);              % y-dist btwn vertices on edge j
 Dzj=(zbj-zaj);              % z-dist btwn vertices on edge j
         
%Calculates intersection point - issue is these are matrices!!
L1=length(Dxi);
L2=length(Dxj);

I1= abs(Dxj.*Dyi-Dxi.*Dyj)>0.0;%1e-15;
I2= (I1==0 & abs(Dyj.*Dzi-Dyi.*Dzj)>0.0);%1e-15
I3= (I1==0 & I2==0);

s=zeros(size(Dxi));
Xstar=zeros(size(Dxi));
Ystar=zeros(size(Dxi));
Zstar=zeros(size(Dxi));
    
if L1==L2

    s(I1)= (Dxj(I1).*(yaj(I1)-yai(I1)) -Dyj(I1).*(xaj(I1)-xai(I1)))./(Dxj(I1).*Dyi(I1)-Dxi(I1).*Dyj(I1)); % problem if s infinite
    
    
    Xstar(I1)=xai(I1)+s(I1).*Dxi(I1);
    Ystar(I1)=yai(I1)+s(I1).*Dyi(I1);
    Zstar(I1)=zai(I1)+s(I1).*Dzi(I1);

    
    
    s(I2)= (Dyj(I2).*(zaj(I2)-zai(I2)) -Dzj(I2).*(yaj(I2)-yai(I2)))./(Dyj(I2).*Dzi(I2)-Dyi(I2).*Dzj(I2)); % problem if s infinite
    Xstar(I2)=xai(I2)+s(I2).*Dxi(I2);
    Ystar(I2)=yai(I2)+s(I2).*Dyi(I2);
    Zstar(I2)=zai(I2)+s(I2).*Dzi(I2);


    
    s(I3)= (Dzj(I3).*(xaj(I3)-xai(I3)) -Dxj(I3).*(zaj(I3)-zai(I3)))./(Dzj(I3).*Dxi(I3)-Dzi(I3).*Dxj(I3)); % problem if s infinite
    Xstar(I3)=xai(I3)+s(I3).*Dxi(I3);
    Ystar(I3)=yai(I3)+s(I3).*Dyi(I3);
    Zstar(I3)=zai(I3)+s(I3).*Dzi(I3);

else
    
    s(I1)= (Dxj.*(yaj-yai(I1)) -Dyj.*(xaj-xai(I1)))./(Dxj.*Dyi(I1)-Dxi(I1).*Dyj); % problem if s infinite
    
    
    Xstar(I1)=xai(I1)+s(I1).*Dxi(I1);
    Ystar(I1)=yai(I1)+s(I1).*Dyi(I1);
    Zstar(I1)=zai(I1)+s(I1).*Dzi(I1);

    
    
    s(I2)= (Dyj.*(zaj-zai(I2)) -Dzj.*(yaj-yai(I2)))./(Dyj.*Dzi(I2)-Dyi(I2).*Dzj); % problem if s infinite
    Xstar(I2)=xai(I2)+s(I2).*Dxi(I2);
    Ystar(I2)=yai(I2)+s(I2).*Dyi(I2);
    Zstar(I2)=zai(I2)+s(I2).*Dzi(I2);


    
    s(I3)= (Dzj.*(xaj-xai(I3)) -Dxj.*(zaj-zai(I3)))./(Dzj.*Dxi(I3)-Dzi(I3).*Dxj); % problem if s infinite
    Xstar(I3)=xai(I3)+s(I3).*Dxi(I3);
    Ystar(I3)=yai(I3)+s(I3).*Dyi(I3);
    Zstar(I3)=zai(I3)+s(I3).*Dzi(I3);
    
end


    
 %Computes Angle using dot product formula to give exterior angle phi and then cos(pi-phi) = cos(Angle) gives interior angle  

 Angle=pi-real(acos((Dxi.*Dxj+Dyi.*Dyj+Dzi.*Dzj)./(sqrt((Dxi.^2)+(Dyi.^2)+(Dzi.^2)).*sqrt((Dxj.^2)+(Dyj.^2)+(Dzj.^2)))));
