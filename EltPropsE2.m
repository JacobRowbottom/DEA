function [a, b, xa, ya, za, xb, yb, zb] = EltPropsE2(xv,yv,zv,sfactx,sfacty,sfactz,CL,step,nSideEls)

		
% Define arclength values at each end of element
% so [ajj,bjj) is a subset of [0,L) where L = CL(nVert)
for j=1:nSideEls       
    a(j,1)=CL+(j-1)*step; 
    b(j,1)=CL+j*step;    
end
% Define the cartesian coords at ajj (or bjj) by moving je-1 (or je)
% steps along the edge jv from the initial vertex in the correct 
% direction of orientation
xa=xv+(a-CL)*sfactx;
xb=xv+(b-CL)*sfactx;
ya=yv+(a-CL)*sfacty;
yb=yv+(b-CL)*sfacty;
za=zv+(a-CL)*sfactz;
zb=zv+(b-CL)*sfactz;

% end of function
end

