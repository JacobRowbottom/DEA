function [xsign ysign slope sfact1 sfact2] = EdgeProps(Lx,Ly)
%Computes some properties of straight edge     


%    Intputs:   
%%%%%    Lx tells you change in x value along edge
%%%%%    Ly tells you change in y value along edge

%    Outputs:   
%%%%%    xsign tells you if edge is increasing (1), decreasing (-1) or constant in the x direction
%%%%%    ysign tells you if edge is increasing (1), decreasing (-1) or constant in the y direction
%%%%%    slope is gradient of edge
%%%%%    sfact1 and sfact2 are factors arising in straight line intersection / parametrisation equations  

if abs(Ly)>0
   
    ysign=Ly/abs(Ly);
else
    
    ysign=0;
end
    
if abs(Lx)>0
		
    slope=Ly/Lx;   %slope of edge jv (when finite)
    sfact1=1/sqrt(1+slope*slope);
    sfact2=abs(slope)*sfact1; 
	xsign=Lx/abs(Lx);
    
else
    
    slope=NaN;
    xsign=0;
	sfact1=0;
	sfact2=1;
    
end

% end of function
end

