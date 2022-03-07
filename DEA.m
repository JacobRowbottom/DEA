%DEA code rewritten using delta direction basis
% Transition from element j to element i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stage=1
eps=1e-15;                          %small parameter near machine precision
epss=1e-8;

%global parameters are abv val of wave vector, damping level and number of vertices of each subdomain

global dAbsP;
global Adamp;
global nVert;  


format long
 params.freq = 50;                                                 % frequency
 params.omega = 2*pi*params.freq;                                    % angular freq = wavenumber if speed=1
 params.rhom = 1/((params.omega)^2);  % density
 params.speed=1;
 params.dAbsP = (1/params.speed);

%om=200*pi; %parameters
%rhof=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometric data for polygonal domains %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nOmega=1;
% % % Define the vertices
% 
% %config A 
% xnodes=(0.25)*[1.4564 0.87 0.0 -0.83 -1.048 -0.28 0.0 0.69];
% ynodes=(0.25)*[0.40381 1.1027 0.6993 1.1720 0.3582 0.0 0.2993 -0.1328];
% %config B
%  xnodes=[1.4564 0.87 0.0 -0.83 -1.048 -0.28 0.0 0.69];
%  ynodes=[0.40381 1.1027 0.9493 1.1720 0.3582 0.0 0.0493 -0.1328];
% %config C
%xnodes=[1.4564 0.87 0 -1.503 -1.503 0 0.69];
%ynodes=[0.40381 1.1027 0.7993 0.7993 0.1993 0.1993 -0.1328];

% % Unit Square
 xnodes=[0 1 1 0];
 ynodes=[0 0 1 1];% 

% % L-shape
% 
% xnodes=[0 1 1 0.5 0.5 0];
% ynodes=[0 0 1 1 0.5 0.5];% 

% xnodes=[0 1 1 1 0.5 0.5 0];
% ynodes=[0 0 0.5 1 1 0.5 0.5];% 


znodes=zeros(size(xnodes));

XYZ=[xnodes' ynodes' znodes'];
 
% Define the subdomains
%

%%Config A and B
% Omega=[5 6 7 3 4; 7 8 1 2 3];

%%Config C
%Omega=[5 6 3 4 0; 6 7 1 2 3];

%%Unit square
Omega=[1 2 3 4];

%%L shape

% Omega=[1 2 5 6; 2 3 4 5];
% Omega = [1 2 3 6 7; 3 4 5 6 0];
%
%for 5-sided polygons
%XYZcent=(XYZ(Omega(:,1),:)+XYZ(Omega(:,2),:)+XYZ(Omega(:,3),:)+XYZ(Omega(:,4),:)+XYZ(Omega(:,5),:))/5;

%for 4 sided polygons
XYZcent=(XYZ(Omega(:,1),:)+XYZ(Omega(:,2),:)+XYZ(Omega(:,3),:)+XYZ(Omega(:,4),:))/4;

%for triangles
%XYZcent=(XYZ(Omega(:,1),:)+XYZ(Omega(:,2),:)+XYZ(Omega(:,3),:))/3;

Lx=zeros(size(Omega));
Ly=zeros(size(Omega));
Lz=zeros(size(Omega));
SL=zeros(size(Omega));
CL=zeros(size(Omega));
Tj=zeros(size(Omega));
Ti=zeros(size(Omega));

%Define global direction set and elements 
%multiples of 4 give symmetry in each quadrant

NQ=2;
N=4*NQ;
dx=1/64;

ThetaC=linspace(0,2*pi,N+1);%endpoints
ThetaAB=0.5*(ThetaC(1:end-1)+ThetaC(2:end));

ThetaStep=ThetaAB(2)-ThetaAB(1);
ThetaC=ThetaC(1:N);

for j=1:nOmega
    nVert(j)= nnz(Omega(j,:));
end

% Vals for splitting into boundary elements
n=zeros(nOmega,1);
nSideEls=zeros(nOmega,max(nVert));
stepb=zeros(nOmega,max(nVert));

tval=0.5; %0.5 is val in my derivation.
WD=1.0; %wall damping factor, value of 1.0 means no losses on reflection

for j=1:nOmega
    %normal vec to each face
    TDNorm(j,:)=cross((XYZ(Omega(j,2),:)-XYZ(Omega(j,1),:)),(XYZ(Omega(j,3),:)-XYZ(Omega(j,1),:)));
    TDNorm(j,:)=TDNorm(j,:)/norm(TDNorm(j,:));

   %vector direction of intersection line between 2 element planes so RotZ and TDNorm coincide
    UVLp=-[TDNorm(j,2) -TDNorm(j,1) 0];
    if norm(UVLp)>=eps
        UVL(j,:)=UVLp/norm(UVLp);
    else
        UVL(j,:)=UVLp;
    end
        
   
    CRA=TDNorm(j,3);
    SRA=sqrt(1-CRA.^2);
    %Defines general rotated x and y axes into plane of element
        RotX(j,:)=[CRA+(1-CRA)*UVL(j,1)^2 UVL(j,1)*UVL(j,2)*(1-CRA)+SRA*UVL(j,3) UVL(j,3)*UVL(j,1)*(1-CRA)-SRA*UVL(j,2)];
        RotY(j,:)=[UVL(j,1)*UVL(j,2)*(1-CRA)-SRA*UVL(j,3) CRA+(1-CRA)*UVL(j,2)^2 UVL(j,3)*UVL(j,2)*(1-CRA)+SRA*UVL(j,1)];
  %      RotZ(j,:)=[UVL(j,1)*UVL(j,3)*(1-CRA)+SRA*UVL(j,2) UVL(j,3)*UVL(j,2)*(1-CRA)-SRA*UVL(j,1) CRA+(1-CRA)*UVL(j,3)^2];

% % %      
    for jj=1:nVert(j)
        jjp=jj+1;
        if jjp>nVert(j)
            jjp=jjp-nVert(j);
        end
               
        %Compute the arc lengths, SL is vector of side lengths and CL is vector of cumulative lengths

        Lx(j,jj)=(xnodes(Omega(j,jjp))-xnodes(Omega(j,jj)));               % x-dist btwn vertices
        Ly(j,jj)=(ynodes(Omega(j,jjp))-ynodes(Omega(j,jj)));               % y-dist btwn vertices
        Lz(j,jj)=(znodes(Omega(j,jjp))-znodes(Omega(j,jj)));               % y-dist btwn vertices
        
        SL(j,jj)=sqrt(((Lx(j,jj))^2)+((Ly(j,jj))^2)+((Lz(j,jj))^2));         % Euclidean dist btwn vertices
        CL(j,jj+1)=CL(j,jj)+SL(j,jj);                         % Cumulative lengths
        
        LEdge(j,jj)=sqrt((Lz(j,jj)*Lz(j,jj))+(Ly(j,jj)*Ly(j,jj))+(Lx(j,jj)*Lx(j,jj)));
      
        sfactx(j,jj)=Lx(j,jj)/LEdge(j,jj);
        sfacty(j,jj)=Ly(j,jj)/LEdge(j,jj);
        sfactz(j,jj)=Lz(j,jj)/LEdge(j,jj);
    end
    
    dAbsP(j)=params.dAbsP;% %0.5*CL(j,nVert(j)+1) Wave-vector 
   % dx= CL(j,nVert(j)+1)/Nelt; %elt size - if want edges discretised same in each subsystem - ignores Nelt and set dx=val.
        %%0.004;%
    %%Helmholtz
    loss=0.01; %hysteretic loss factor
    Adamp(j)=loss*params.omega*dAbsP(j)/2;
    
    for iv=1:nVert(j)
        nSideEls(j,iv)=max(1,round(SL(j,iv)/dx)); %No. of elements on each edge
        stepb(j,iv)=SL(j,iv)/nSideEls(j,iv); %step size on each edge     
        n(j)=n(j)+nSideEls(j,iv);              %total number of boundary elements
    end
end


OmegaE=zeros(size(Omega));% gives global edge number associated to the pair (j,jv)
EdgeOK=(zeros(sum(nVert)));
ACheck=(zeros(nOmega,max(nVert)));
for j=1:nOmega  
    % defines self-element interactions
   EdgeOK(1+sum(nVert(1:j-1)):nVert(j)+sum(nVert(1:j-1)),1+sum(nVert(1:j-1)):nVert(j)+sum(nVert(1:j-1)))=ones(nVert(j))-eye(nVert(j)); 

   for jv=1:nVert(j)     
       jvp=jv+1;     
       while jvp>nVert(j)
             jvp=jvp-nVert(j);
       end
       jvpp=jv+2;
       while jvpp>nVert(j)
             jvpp=jvpp-nVert(j);
       end
       
       PListP(:,jv,j)=[Omega(j,jv); Omega(j,jvp)];
       %changed for mesh case - must change back in 2D general domains
       % ACheck(j,jv)=det([1 xnodes(Omega(j,jv)) ynodes(Omega(j,jv));1  xnodes(Omega(j,jvp)) ynodes(Omega(j,jvp));1  xnodes(Omega(j,jvpp)) ynodes(Omega(j,jvpp)) ]);
       ACheck(j,jv)=1; % simplification for 3D surface mesh 
       
       if ACheck(j,jv)==0 % removes collinear mapping
          
           EdgeOK(jv+sum(nVert(1:j-1)),jvp+sum(nVert(1:j-1)))=0;
           EdgeOK(jvp+sum(nVert(1:j-1)),jv+sum(nVert(1:j-1)))=0;
           
       end
%       for jE=1:length(Edges)
%           if (Edges(jE,1)==Omega(j,jv)  && Edges(jE,2)==Omega(j,jvp)) || (Edges(jE,1)==Omega(j,jvp)  && Edges(jE,2)==Omega(j,jv))
%              OmegaE(j,jv)=jE;
%           end      
%       end       
   end
end
FreeEdges=ones(nOmega,max(nVert));

for j=1:nOmega
    j
    for jv=1:nVert(j)
        jvp=jv+1;     
        while jvp>nVert(j)
              jvp=jvp-nVert(j);
        end
        jvm=jv-1;     
        while jvm<1
              jvm=jvm+nVert(j);
        end
        
        for i=1:nOmega
            if j~=i
                for iv=1:nVert(i)
                    if PListP(1:2,jv,j)==PListP(2:-1:1,iv,i) % tests if matching edges for mapping to connected edges                       
                       
                        EdgeOK(iv+sum(nVert(1:i-1)),sum(nVert(1:j-1))+1:sum(nVert(1:j-1))+nVert(j))=ones(1,nVert(j));               
                        EdgeOK(iv+sum(nVert(1:i-1)),jv+sum(nVert(1:j-1)))=0;                    
                        if ACheck(j,jv)==0           %prevents colinear transmission            
                            EdgeOK(iv+sum(nVert(1:i-1)),jvp+sum(nVert(1:j-1)))=0;                             
                        end
                        if ACheck(j,jvm)==0                       
                            EdgeOK(iv+sum(nVert(1:i-1)),jvm+sum(nVert(1:j-1)))=0;                             
                        end                    
                        FreeEdges(j,jv)=0; 
                    end
                end
            end
        end
    end
end

%% source Vector data
Source=0; %0 for bc, 1 for point sc
ScDomVert=[];
%ScVert=1;
if Source==0
    % Case 1: BC along edge x=a;
    a=min(xnodes);
    strength=1;%params.omega;
    t0=(pi/4);% %source direction
%Helmholtz
%     PreFactB=cos(t0)*params.rhom/((1+(loss^2)/16));
    %Biharmonic
    PreFactB=(strength)*cos(t0);
    else
    % Case 2: pt source at (x,y)=(Ox,Oy);
    ScDom=1; % pick a sub-domain for the source point
%Conf C  
% Ox=-1.4;
% Oy=0.4993;
%Conf A and B
%  Ox=-0.1;
%  Oy=0.125;
 %unit square centre
Ox=0.5;
Oy=0.5;
% Other positions
%  Ox=0.25;
%  Oy=0.25;
 
  Oz=0.0;
    %%Helmholtz
    PreFact=params.rhom*params.omega*params.omega/(8*pi*params.omega);%* - include these for energy density params.rhom*params.omega*params.omega
end

%%Final density preprocessing:
  addpath(genpath('distmesh'));
  PpOm=0.025; 
% 
stage=2
%% Starts main loop over edges (equiv vertices) jv transmitting from 
%Convention: an edge is labelled by its 1st vertex with anti-clockwise orientation
for j=1:nOmega   
    
    %preprocess for final density
 % if general polygons  
       XN=xnodes(Omega(j,1:nVert(j)))';
       YN=ynodes(Omega(j,1:nVert(j)))';
       pv=[XN YN;XN(1) YN(1)]; 
% % %     figure(j)
% % %     
      [points,tri]=distmesh2d(@dpoly,@huniform,PpOm,[min(XN),min(YN); max(XN),max(YN)],pv,pv);
% % %   
       XP=(points(tri(:,1),1)+points(tri(:,2),1)+points(tri(:,3),1))/3;
       YP=(points(tri(:,1),2)+points(tri(:,2),2)+points(tri(:,3),2))/3;
       ZP=zeros(size(XP));
     
    PlotSize(j)=length(XP);
    DirSplit=zeros(max(nVert),PlotSize(j));
    
    FDValV=[];%zeros(2*NQ,max(nVert));
    FDPosV=[];%zeros(2*NQ,max(nVert));
    LGStoreV=[];%zeros(2*NQ,max(nVert));
    FDScaleV=[];
   
    for pts=1:PlotSize(j)
        for jv=1:nVert(j)
% %         % Gives end-point directions to split global directions into those hitting each individual edge of
% %         % the subdomain j from (XP,YP,ZP)
         
            DirVec=XYZ(Omega(j,jv),:)-[XP(pts) YP(pts) ZP(pts)];
            DirX=dot(DirVec,RotX(j,:));
            DirY=dot(DirVec,RotY(j,:));
            
            DirSplit(jv,pts)=atan2(DirY,DirX);
            while DirSplit(jv,pts)<0
                  DirSplit(jv,pts)=2*pi+DirSplit(jv,pts);
            end
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for jv=1:nVert(j)    
        
        jvp=jv+1;
         
         while jvp>nVert(j)
            jvp=jvp-nVert(j);
         end
                                   
         xva=xnodes(Omega(j,jv));
         xvb=xnodes(Omega(j,jvp));
         yva=ynodes(Omega(j,jv));
         yvb=ynodes(Omega(j,jvp));  
         zva=znodes(Omega(j,jv));
         zvb=znodes(Omega(j,jvp));
         
         EdgeX=dot([xvb-xva yvb-yva zvb-zva],RotX(j,:));
         EdgeY=dot([xvb-xva yvb-yva zvb-zva],RotY(j,:));
        % EdgeZ=dot([xvb-xva yvb-yva zvb-zva],RotZ(j,:));
         
         lowercut=atan2(EdgeY,EdgeX);
        
         while lowercut<0
               lowercut=lowercut+2*pi; 
         end
        
        uppercut=lowercut+pi;
        
        while uppercut>=2*pi
              uppercut=uppercut-2*pi;
        end
        
        if lowercut<uppercut      
            ThetaI=find(ThetaC>=lowercut & ThetaC<uppercut);
            ThetaLoc=0.5*pi-ThetaC(ThetaI)+lowercut;
            
        else     
            ThetaI=[find(ThetaC>=lowercut) find(ThetaC<uppercut)];
            ThetaLocL=0.5*pi-ThetaC(ThetaI)+lowercut;
            ThetaLocU=-0.5*pi-ThetaC(ThetaI)+uppercut;
            ThetaLocL(abs(ThetaLocL)>0.5*pi)=[];
            ThetaLocU(abs(ThetaLocU)>0.5*pi)=[];
            ThetaLoc=[ThetaLocL ThetaLocU];
           
        end

        Tj(j,jv)=length(ThetaI);
      
        alpha=atan2(EdgeX,-EdgeY);
      
        %slopeij=tan(alpha-ThetaLoc);
        %ThetaLoc
        
        %SlopeMat=repmat(slopeij',[1,nSideEls(j,jv)]);              
        %GlobMat=ThetaC(ThetaI)';%repmat(,[1,nSideEls(j,jv)]);
        
        %blocks correspond to number of local directions for each element
        %on edge jv: size is 2NQ directions per element
        jBlockA=(sum(sum(Tj(1:j-1,1:end).*nSideEls(1:j-1,1:end)))+sum(Tj(j,1:jv-1).*nSideEls(j,1:jv-1)))+1;
       
        jBlockB=(sum(sum(Tj(1:j-1,1:end).*nSideEls(1:j-1,1:end)))+sum(Tj(j,1:jv).*nSideEls(j,1:jv)));
        
        %DSj=stepb(j,jv)*ones(1,Tj(j,jv)*nSideEls(j,jv));
               
        %DSjMat=reshape(DSj,[Tj(j,jv),nSideEls(j,jv)]);
              
       
        [ajj, bjj, xaj, yaj, zaj, xbj, ybj, zbj] = EltPropsE2(xva,yva,zva,sfactx(j,jv),sfacty(j,jv),sfactz(j,jv),CL(j,jv),stepb(j,jv),nSideEls(j,jv));
        
        % for the interior density %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        PlotCount=0;%zeros(PlotSize(j),1);
       
      
         for pts=1:PlotSize(j)  
             %FDVal=zeros(Tj(j,jv),1);
             ScaleFact=zeros(Tj(j,jv),1);
             FDPos=zeros(Tj(j,jv),1);
             LGStore=zeros(Tj(j,jv),1);
             
             sElt=zeros(nSideEls(j,jv),Tj(j,jv));
             tList1=zeros(Tj(j,jv),Tj(j,jv));
             sEltList=zeros(Tj(j,jv),Tj(j,jv));
             LGam=zeros(Tj(j,jv),1);
             xGam=0.*ones(Tj(j,jv),1);
             yGam=0.*ones(Tj(j,jv),1);
             zGam=0.*ones(Tj(j,jv),1);
             AngleGam=zeros(Tj(j,jv),1);
             sGam=zeros(Tj(j,jv),1);
             ThetaGam=zeros(Tj(j,jv),1);

            
             if DirSplit(jvp,pts)>=DirSplit(jv,pts) %use DirSplit to separate thetac into directions for given edge jv to jvp
                 
                 DirEdge=ThetaC(DirSplit(jv,pts)<=ThetaC & ThetaC<=DirSplit(jvp,pts));%+2*eps-2*eps
             else
                 DirEdge=ThetaC(DirSplit(jv,pts)<=ThetaC | ThetaC<=DirSplit(jvp,pts));%+2*eps-2*eps
             end             
           
            
             Ix=XP(pts)+cos(DirEdge)*RotX(j,1)+sin(DirEdge)*RotY(j,1); %Another point on the ray moving along a unit vector in direction DirEdge
             Iy=YP(pts)+cos(DirEdge)*RotX(j,2)+sin(DirEdge)*RotY(j,2);
             Iz=ZP(pts)+cos(DirEdge)*RotX(j,3)+sin(DirEdge)*RotY(j,3);
             %(xGam, yGam) are points on gamma reached by all rays from solution point
             %     
             
                        
             for kk=1:length(Ix)
               
                 [xGam(kk), yGam(kk), zGam(kk), AngleGam(kk)] = Intersection3D(Ix(kk),Iy(kk),Iz(kk),XP(pts),YP(pts),ZP(pts),xva,yva,zva,xvb,yvb,zvb);
                 
                 sGam(kk)=CL(j,jv)+sqrt((xGam(kk)-xva).^2+(yGam(kk)-yva).^2+(zGam(kk)-zva).^2); %distance of each solution point trajectory along edge jv
                 
                 ThetaGam(kk)=AngleGam(kk)-0.5*pi;
               %  pause
                 
                 LGam(kk)=sqrt((xGam(kk)-XP(pts)).^2 + (yGam(kk)-YP(pts)).^2+(zGam(kk)-ZP(pts)).^2);
                 
                               
                 sElt(:,kk)=bsxfun(@ge,sGam(kk)*ones(nSideEls(j,jv),1),ajj); %checks if sGam greater than elt starting with ajj
                                  
                 tList1(1:Tj(j,jv),kk)=bsxfun(@lt,abs(ThetaGam(kk)*ones(Tj(j,jv),1)-ThetaLoc'),10*eps); %10*eps check if thetagam in local set of edge angles
                     
                 sEltList(:,kk)=sum(sElt(:,kk),'double'); %1 column for each edge element. Sum includes a 1 for each element sGam is greater than.
                 %Need to think carefully how this bit works!
             end
          
             KeyFD=tList1.*sEltList; %final density
              
             [KRow, KCol, KVal]=find(KeyFD);
             
            
             FDPos(1:length(KVal),1)=(jBlockA-1)*ones(size(KVal))+Tj(j,jv)*(KVal-ones(size(KVal)))+KRow;
           
             

            % FDVal(1:length(KVal),1)=DirEdge(KCol); 
             LGStore(1:length(KVal),1)=LGam(KCol); 
             
            % ScaleV(((jv-1)*2*NQ+1:jv*2*NQ))=DSjMat(:,KVal);
             FDPosV((jv-1)*Tj(j,jv)+1:jv*Tj(j,jv),pts)=FDPos;
            % FDValV((jv-1)*Tj(j,jv)+1:jv*Tj(j,jv),pts)=FDVal;
            
            FDScaleV((jv-1)*Tj(j,jv)+1:jv*Tj(j,jv),pts)=(stepb(j,jv)^(-tval))./cos(ThetaGam);
%           

                   
             LGStoreV((jv-1)*Tj(j,jv)+1:jv*Tj(j,jv),pts)=LGStore;

         end
        
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        %%%%%%%%%%%SOURCE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if Source==0         
            % case 1: edge BC along any line with 0 x-coordinate along its length 
            ScVec(jBlockA:jBlockB,1)=zeros(Tj(j,jv)*nSideEls(j,jv),1);
            midpart=zeros(nSideEls(j,jv),1);
            ScDir=zeros(Tj(j,jv),1);

            if abs(xva-a)<eps && abs(xvb-a)<eps % && abs(ybj)<(0.75+eps) && abs(yaj)>(0.25-eps)
           % if j==1 && jv==4
                % for direction =t0 rads.  
            
                ScDir(abs(ThetaC(ThetaI)-t0)<=0.51*ThetaStep)=1;
                ScDir(abs(ThetaC(ThetaI)-(t0+2*pi))<=0.51*ThetaStep)=1;% in case t0=0;
                 ScDir=repmat(ScDir,nSideEls(j,jv),1);
%                 midpart(abs(yaj)<=(0.726833761643633) & abs(ybj)>(0.2694-eps))=1; %% nb - yaj>ybj!
%                 
%                ScDir=kron(midpart,ScDir);
                  
                ScVec(jBlockA:jBlockB,1)=sparse((PreFactB)*ScDir.*(stepb(j,jv)^tval));   %  /SL(j,jv)  
            end
            
        else
            %%
        % Case 2: pt source at (x,y)=(Ox,Oy,Oz);
        
            ScInt=zeros(Tj(j,jv),nSideEls(j,jv));
            
            
            % first consider reflected contributions
            % 2nd and 3rd conditions below not needed for source inside domain
            if any(j==ScDom) %&& ScVert~=PListP(1,jv,j) && ScVert~=PListP(2,jv,j)% restrict to source point domains
            %need to integrate ScTerm entries and place into right positions in ScVec
            
               if FreeEdges(j,jv)==0
                        
                   for loop1=1:Tj(j,jv)
                       for loop2=1:nSideEls(j,jv)                  
                                               
                            SR = @(s) SourceIntegrandR(s,xaj(loop2),yaj(loop2),zaj(loop2),xbj(loop2),ybj(loop2),zbj(loop2),Ox,Oy,Oz,Adamp(j),ThetaLoc(loop1),ThetaStep);
         
                            ScInt(loop1,loop2)=integral(SR,0,stepb(j,jv));
                            
                                                     
                       end
                   end
               else
                        
                  for loop1=1:Tj(j,jv)
                     for loop2=1:nSideEls(j,jv)                  
                                      %reflecting outer edge         
                          SE = @(s) SourceIntegrandE(s,xaj(loop2),yaj(loop2),zaj(loop2),xbj(loop2),ybj(loop2),zbj(loop2),Ox,Oy,Oz,Adamp(j),ThetaLoc(loop1),ThetaStep);
                                
            %WD = wall damping factor
                          ScInt(loop1,loop2)=WD*integral(SE,0,stepb(j,jv));  %   ,'abstol',1e-12,'reltol',1e-8                                
                     
                     end
                  end      
               end
             
            end           
            % now consider transmitted contributions
             for sv=1:nVert(ScDom) % for source inside domain
            % for sv=1:length(ScDom)  % for vertex source
               
                 if PListP(1:2,sv,ScDom)==PListP(2:-1:1,jv,j) %for source in domain:  restrict to source domain connecting edges 
              %    if PListP(1:2,ScDomVert(sv,2),ScDom(sv))==PListP(2:-1:1,jv,j)  %for source on vertex
               %                      
               
                    for loop1=1:Tj(j,jv)
                         for loop2=1:nSideEls(j,jv)                  
%                           
%source on vertex:
     %                        ST = @(s) SourceIntegrandT(s,xaj(loop2),yaj(loop2),zaj(loop2),xbj(loop2),ybj(loop2),zbj(loop2),Ox,Oy,Oz,Adamp(ScDom(sv)),ThetaLoc(loop1),ThetaStep);
  %source in subdomain:    
                        ST = @(s) SourceIntegrandT(s,xaj(loop2),yaj(loop2),zaj(loop2),xbj(loop2),ybj(loop2),zbj(loop2),Ox,Oy,Oz,Adamp(ScDom),ThetaLoc(loop1),ThetaStep);

                             %                          
                             ScInt(loop1,loop2)=integral(ST,0,stepb(j,jv));         
                         end
                    end
                 end
              end
%             if nnz(ScInt(:))>0
%                 ScInt
%                 pause
%             end
            ScVec(jBlockA:jBlockB,1)=sparse(PreFact.*ScInt(:));
            
        end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
           
        for i=1:nOmega   
            for iv=1:nVert(i)    
                          
                   ivp=iv+1;
         
                    while ivp>nVert(i)
                        ivp=ivp-nVert(i);
                    end
                                
                    xvai=xnodes(Omega(i,iv));
                    xvbi=xnodes(Omega(i,ivp));
                    yvai=ynodes(Omega(i,iv));
                    yvbi=ynodes(Omega(i,ivp));
                    zvai=znodes(Omega(i,iv));
                    zvbi=znodes(Omega(i,ivp));
         
                    TurnAxis=[xvbi-xvai yvbi-yvai zvbi-zvai];
                    
                    EdgeXi=dot(TurnAxis,RotX(i,:));
                    EdgeYi=dot(TurnAxis,RotY(i,:));
                    % EdgeZ=dot([xvb-xva yvb-yva zvb-zva],RotZ(j,:));
         
                    lowercuti=atan2(EdgeYi,EdgeXi);
                            
                    while lowercuti<0
                        lowercuti=lowercuti+2*pi; 
                    end
            
                    uppercuti=lowercuti+pi;
        
                    while uppercuti>=2*pi
                        uppercuti=uppercuti-2*pi;
                    end
        
                    if lowercuti<uppercuti      
                        ThetaIi=find(ThetaC>=lowercuti & ThetaC<uppercuti);
                        %ThetaLoci=0.5*pi-ThetaC(ThetaIi)+lowercuti;
            
                    else     
                        ThetaIi=[find(ThetaC>=lowercuti) find(ThetaC<uppercuti)];
                       % ThetaLocLi=0.5*pi-ThetaC(ThetaIi)+lowercuti;
                        %ThetaLocUi=-0.5*pi-ThetaC(ThetaIi)+uppercuti;
                       % ThetaLocLi(abs(ThetaLocLi)>0.5*pi)=[];
                      %  ThetaLocUi(abs(ThetaLocUi)>0.5*pi)=[];
                       % ThetaLoci=[ThetaLocLi ThetaLocUi];
                    end
                    Ti(i,iv)=length(ThetaIi);
                    alphai=atan2(EdgeXi,-EdgeYi);
                    
                    %AlphaIMat=alphai*ones(Ti(i,iv),Tj(j,jv));
                    %AlphaMat=alpha*ones(Ti(i,iv),Tj(j,jv));
                    ThetaMat=kron(ones(nSideEls(j,jv),1),ThetaC(ThetaI));
                    ThetaIMat=kron(ones(nSideEls(i,iv),1),ThetaC(ThetaI));
                   
                    
                    iBlockA=(sum(sum(Ti(1:i-1,1:end).*nSideEls(1:i-1,1:end)))+sum(Ti(i,1:iv-1).*nSideEls(i,1:iv-1)))+1;
                    iBlockB=(sum(sum(Ti(1:i-1,1:end).*nSideEls(1:i-1,1:end)))+sum(Ti(i,1:iv).*nSideEls(i,1:iv)));
                    
                 if EdgeOK(iv+sum(nVert(1:i-1)),jv+sum(nVert(1:j-1)))==1 %avoids parallel or non-connecting edges
                   % DSi=stepb(i,iv)*ones(Ti*nSideEls(i,iv),1);
                    
                   %%update
                   RT=ones(Ti(i,iv),Tj(j,jv));
                   %% 
                   
                   %DSiMat=reshape(DSi,[Ti(i,iv),nSideEls(i,iv)]);
       
                    StepProd=(stepb(i,iv)^(-0.5))*(stepb(j,jv)^(-0.5));%.*SL(i,iv)/SL(j,jv);
                                                               
                    Bp=zeros(nSideEls(i,iv),Tj(j,jv)*nSideEls(j,jv));
                 %   slopeiji=tan(alphai-ThetaLoci);

                    [~,~, xai, yai, zai, xbi, ybi, zbi] = EltPropsE2(xvai,yvai,zvai,sfactx(i,iv),sfacty(i,iv),sfactz(i,iv),CL(i,iv),stepb(i,iv),nSideEls(i,iv));
          
  %%%%%%%%%%%%%%%    NEED TO CORRECT FROM HERE !!
                                     
     %               SlopeMatV=SlopeMat(:);
                   % GlobMatV=GlobMat(:);

%                      xaiMat=kron(xai,ones(1,nSideEls(j,jv)));
%                      yaiMat=kron(yai,ones(1,nSideEls(j,jv)));
%                      zaiMat=kron(zai,ones(1,nSideEls(j,jv)));
%                      xbiMat=kron(xbi,ones(1,nSideEls(j,jv)));
%                      ybiMat=kron(ybi,ones(1,nSideEls(j,jv)));
%                      zbiMat=kron(zbi,ones(1,nSideEls(j,jv)));
% %                     
                      xajMat=(xaj*ones(1,Tj(j,jv)))';
                      yajMat=(yaj*ones(1,Tj(j,jv)))';
                      zajMat=(zaj*ones(1,Tj(j,jv)))';
% % %                    
                       xbjMat=(xbj*ones(1,Tj(j,jv)))';
                       ybjMat=(ybj*ones(1,Tj(j,jv)))';
                       zbjMat=(zbj*ones(1,Tj(j,jv)))';
%                     
      %              SlopeMatK=kron(ones(2*NQ*nSideEls(i,iv),1),SlopeMatV');
                   
                   
                    %Another point on ray backprojected from (xai,yai,zai)
                    XRA=xai+cos(ThetaIMat)*RotX(j,1)+sin(ThetaIMat)*RotY(j,1); %Another point on the ray moving along a unit vector in direction DirEdge
                    YRA=yai+cos(ThetaIMat)*RotX(j,2)+sin(ThetaIMat)*RotY(j,2);
                    ZRA=zai+cos(ThetaIMat)*RotX(j,3)+sin(ThetaIMat)*RotY(j,3);

                    %Another point on ray backprojected from (xbi,ybi,zbi)
                    XRB=xbi+cos(ThetaIMat)*RotX(j,1)+sin(ThetaIMat)*RotY(j,1); %Another point on the ray moving along a unit vector in direction DirEdge
                    YRB=ybi+cos(ThetaIMat)*RotX(j,2)+sin(ThetaIMat)*RotY(j,2);
                    ZRB=zbi+cos(ThetaIMat)*RotX(j,3)+sin(ThetaIMat)*RotY(j,3);
                   
                    
                    Xa=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    Ya=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    Za=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    Xb=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    Yb=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    Zb=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    DaAj=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    DbAj=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    DaBj=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    DbBj=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    DaCj=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    DbCj=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                   
                    xcj=0.5*(xaj+xbj);
                    ycj=0.5*(yaj+ybj);
                    zcj=0.5*(zaj+zbj);
                    
                    for ie=1:nSideEls(j,jv)
                        
                        [Xa(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv)), Ya(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv)), Za(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv)), ~] = Intersection3D(XRA,YRA,ZRA,xai,yai,zai,xaj(ie),yaj(ie),zaj(ie),xbj(ie),ybj(ie),zbj(ie));
                        [Xb(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv)), Yb(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv)), Zb(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv)), ~] = Intersection3D(XRB,YRB,ZRB,xbi,ybi,zbi,xaj(ie),yaj(ie),zaj(ie),xbj(ie),ybj(ie),zbj(ie));
                   
                        DaAj(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))=sqrt((Xa(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-xaj(ie)).^2+(Ya(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-yaj(ie)).^2+(Za(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-zaj(ie)).^2);
                        DbAj(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))=sqrt((Xb(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-xaj(ie)).^2+(Yb(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-yaj(ie)).^2+(Zb(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-zaj(ie)).^2);
                        DaBj(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))=sqrt((Xa(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-xbj(ie)).^2+(Ya(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-ybj(ie)).^2+(Za(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-zbj(ie)).^2);
                        DbBj(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))=sqrt((Xb(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-xbj(ie)).^2+(Yb(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-ybj(ie)).^2+(Zb(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-zbj(ie)).^2);
                       
                        DaCj(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))=sqrt((Xa(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-xcj(ie)).^2+(Ya(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-ycj(ie)).^2+(Za(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-zcj(ie)).^2);
                        DbCj(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))=sqrt((Xb(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-xcj(ie)).^2+(Yb(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-ycj(ie)).^2+(Zb(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-zcj(ie)).^2);                 
                    
                    end
                    
                    AA=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    BB=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    BpU=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    BpL=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    
                    %????????????????
                         XAj=repmat(xajMat(:)',nSideEls(i,iv),1);
                         YAj=repmat(yajMat(:)',nSideEls(i,iv),1);
                         ZAj=repmat(zajMat(:)',nSideEls(i,iv),1);
                         XBj=repmat(xbjMat(:)',nSideEls(i,iv),1);
                         YBj=repmat(ybjMat(:)',nSideEls(i,iv),1);
                         ZBj=repmat(zbjMat(:)',nSideEls(i,iv),1);
                        
                                       
                    AA(DaAj>DaBj)=DaCj(DaAj>DaBj);
                    AA(DaAj<DaBj)=-DaCj(DaAj<DaBj);
                    
                    BB(DbAj>DbBj)=DbCj(DbAj>DbBj);
                    BB(DbAj<DbBj)=-DbCj(DbAj<DbBj);

                  
                    BpU(AA>BB)=bsxfun(@min,0.5*stepb(j,jv),AA(AA>BB));
                    BpL(AA>BB)=bsxfun(@max,-0.5*stepb(j,jv),BB(AA>BB));                       
                    BpU(BB>AA)=bsxfun(@min,0.5*stepb(j,jv),BB(BB>AA));
                    BpL(BB>AA)=bsxfun(@max,-0.5*stepb(j,jv),AA(BB>AA));
                    
                    %AA,BB correspond to the back-projection of the ray and 
                    %corresponds to integration limit inside the element,
                    %+ or - 0.5DSj is + or - half element length and
                    %coresponds to (xaj,yaj).
                 
                    if Adamp(j)<eps
                        
                                            
                        Bp(BpU>BpL)=BpU(BpU>BpL)-BpL(BpU>BpL);
                            
                    else
                                               
                        %back proj orientations preserve direction of
                        %destination edge
                        
                         % when preimage outside element, set XBj, YBj to be
                        % lower element bdry and A to be upper
                        
                        
                        if i==j
                             XAj(BpL<BpU & abs(BpL+0.5*stepb(j,jv))>=eps)=Xb(BpL<BpU & abs(BpL+0.5*stepb(j,jv))>=eps);
                             YAj(BpL<BpU & abs(BpL+0.5*stepb(j,jv))>=eps)=Yb(BpL<BpU & abs(BpL+0.5*stepb(j,jv))>=eps);                           
                             ZAj(BpL<BpU & abs(BpL+0.5*stepb(j,jv))>=eps)=Zb(BpL<BpU & abs(BpL+0.5*stepb(j,jv))>=eps);                           

                             XBj(BpL<BpU & abs(BpU-0.5*stepb(j,jv))>=eps)=Xa(BpL<BpU & abs(BpU-0.5*stepb(j,jv))>=eps);
                             YBj(BpL<BpU & abs(BpU-0.5*stepb(j,jv))>=eps)=Ya(BpL<BpU & abs(BpU-0.5*stepb(j,jv))>=eps);                   
                             ZBj(BpL<BpU & abs(BpU-0.5*stepb(j,jv))>=eps)=Za(BpL<BpU & abs(BpU-0.5*stepb(j,jv))>=eps);                   
                                                                       
                        else
                             XAj(BpL<BpU & abs(BpL+0.5*stepb(j,jv))>=eps)=Xa(BpL<BpU & abs(BpL+0.5*stepb(j,jv))>=eps);
                             YAj(BpL<BpU & abs(BpL+0.5*stepb(j,jv))>=eps)=Ya(BpL<BpU & abs(BpL+0.5*stepb(j,jv))>=eps);                           
                             ZAj(BpL<BpU & abs(BpL+0.5*stepb(j,jv))>=eps)=Za(BpL<BpU & abs(BpL+0.5*stepb(j,jv))>=eps);                           

                             XBj(BpL<BpU & abs(BpU-0.5*stepb(j,jv))>=eps)=Xb(BpL<BpU & abs(BpU-0.5*stepb(j,jv))>=eps);
                             YBj(BpL<BpU & abs(BpU-0.5*stepb(j,jv))>=eps)=Yb(BpL<BpU & abs(BpU-0.5*stepb(j,jv))>=eps);                   
                             ZBj(BpL<BpU & abs(BpU-0.5*stepb(j,jv))>=eps)=Zb(BpL<BpU & abs(BpU-0.5*stepb(j,jv))>=eps);                   
                        end
                   GlobMat=repmat(ThetaIMat,1,nSideEls(j,jv));
 
                        %Another point on ray backprojected from (xai,yai,zai)
                    XRAi=XAj+cos(GlobMat)*RotX(j,1)+sin(GlobMat)*RotY(j,1); %Another point on the ray moving along a unit vector in direction DirEdge
                    YRAi=YAj+cos(GlobMat)*RotX(j,2)+sin(GlobMat)*RotY(j,2);
                    ZRAi=ZAj+cos(GlobMat)*RotX(j,3)+sin(GlobMat)*RotY(j,3);
                    
                    %Another point on ray backprojected from (xbi,ybi,zbi)
                    XRBi=XBj+cos(GlobMat)*RotX(j,1)+sin(GlobMat)*RotY(j,1); %Another point on the ray moving along a unit vector in direction DirEdge
                    YRBi=YBj+cos(GlobMat)*RotX(j,2)+sin(GlobMat)*RotY(j,2);
                    ZRBi=ZBj+cos(GlobMat)*RotX(j,3)+sin(GlobMat)*RotY(j,3);
                    
                    XAi=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    YAi=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    ZAi=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    XBi=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    YBi=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    ZBi=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));

                    for ie=1:nSideEls(i,iv)
                        [XAi(ie,:), YAi(ie,:), ZAi(ie,:), dummy] = Intersection3D(XRAi(ie,:),YRAi(ie,:),ZRAi(ie,:),XAj(ie,:),YAj(ie,:),ZAj(ie,:),xai(ie),yai(ie),zai(ie),xbi(ie),ybi(ie),zbi(ie));
                        [XBi(ie,:), YBi(ie,:), ZBi(ie,:), dummy] = Intersection3D(XRBi(ie,:),YRBi(ie,:),ZRBi(ie,:),XBj(ie,:),YBj(ie,:),ZBj(ie,:),xai(ie),yai(ie),zai(ie),xbi(ie),ybi(ie),zbi(ie));
                    end
                        % for paralel edges will have problem with La=Lb
                       
                        
                        Lb=sqrt((XBi-XBj).^2+(YBi-YBj).^2+(ZBi-ZBj).^2);
                        La=sqrt((XAi-XAj).^2+(YAi-YAj).^2+(ZAi-ZAj).^2);
                        dLds=(Lb-La)./(BpU-BpL);
             
       
                      % plot(Lj(BpU>BpL)-BpU(BpU>BpL)+BpL(BpU>BpL))
                       
                       % pause
                       
                        if La(BpU>BpL)~=Lb(BpU>BpL)
                                                        
                            Bp(BpU>BpL)=(-exp(-Adamp(j)*Lb(BpU>BpL))+exp(-Adamp(j)*La(BpU>BpL)))./(dLds(BpU>BpL)*Adamp(j));                       
                                                        
                        else %parallel edge case, L is const
                            
                            Bp(BpU>BpL)=(exp(-Adamp(j)*Lb(BpU>BpL))).*(BpU(BpU>BpL)-BpL(BpU>BpL));

                        end
                       
                              
                    end
                  
                    ThetaCMat=kron(ThetaC(ThetaIi)',ones(1,Tj(j,jv)));
                  
                    %ThetaIiMat=kron(ThetaIi',ones(1,2*NQ));
                    %ThetaIMat=kron(ones(2*NQ,1),ThetaI);
                         
                    DirMap=zeros(Ti(i,iv),Tj(j,jv));
                   
                     if i==j

                        AngleRefIn=pi+ThetaC(ThetaI);
                        AngleRefOut=mod(2*alphai-AngleRefIn,2*pi);          
%                         pause
                        AngleRefM=kron(ones(Ti(i,iv),1),AngleRefOut);
                                          
                                             
                        if FreeEdges(i,iv)==0

                          
                            RT=0*ones(Ti(i,iv),Tj(j,jv));     
                               DirMap(abs(AngleRefM-ThetaCMat)<0.5*ThetaStep)=1;%
                           
                                             
                         else             %WD wall damping factor          
                           RT= WD*ones(Ti(i,iv),Tj(j,jv));
                                       
                           DirMap(abs(AngleRefM-ThetaCMat)<0.5*ThetaStep)=1;%
                           
                        end
                        
                        DirMap(mod(abs(ThetaCMat-alphai+0.5*pi),2*pi)==0)=0;
                        ThetaIIMat=kron(ones(Ti(i,iv),1),ThetaC(ThetaI));
                        DirMap(mod(abs(ThetaIIMat-alpha+0.5*pi),2*pi)==0)=0;

                      
                    else
                                          
                        TurnAng=acos(dot(TDNorm(i,:),TDNorm(j,:)));
                                  
                        TurnAxisU=-TurnAxis./norm(TurnAxis);
                     
                        ConTest=dot(XYZcent(i,:)-XYZcent(j,:),TDNorm(i,:)-TDNorm(j,:));
                        %Tests if faces convex of concave so rotate about
                        %correct angle
                        
                        if ConTest<0 % correction for concave edges
                            TurnAng=-TurnAng;
                        end
                        
                        cT=cos(TurnAng);
                        sT=sin(TurnAng);
                        %TurnX is the local axis in element j, rotated into the plane of element i   
                        TurnX=[cT+(1-cT)*TurnAxisU(1)^2 (1-cT)*TurnAxisU(1)*TurnAxisU(2)-sT*TurnAxisU(3) (1-cT)*TurnAxisU(1)*TurnAxisU(3)+sT*TurnAxisU(2); (1-cT)*TurnAxisU(1)*TurnAxisU(2)+sT*TurnAxisU(3) cT+(1-cT)*TurnAxisU(2)^2 (1-cT)*TurnAxisU(2)*TurnAxisU(3)-sT*TurnAxisU(1); (1-cT)*TurnAxisU(1)*TurnAxisU(3)-sT*TurnAxisU(2) (1-cT)*TurnAxisU(2)*TurnAxisU(3)+sT*TurnAxisU(1) cT+(1-cT)*TurnAxisU(3)^2]*(RotX(j,:)');
                        
                        SignTest=dot(TurnX,RotY(i,:));

                        % Rotates the directions in element j to orient with the
                        % local x-axis in element i
                        AngleTran=mod(ThetaC(ThetaI)+sign(SignTest)*abs(acos(min(dot(RotX(i,:),TurnX),1))),2*pi);%
                         
                        AngleTran(abs(AngleTran-2*pi)<2e-8)=0;%AngleTran(abs(AngleTran-2*pi)<2e-8)-AngleTran(abs(AngleTran-2*pi)<2e-8);

                        AngleTranM=kron(ones(Ti(i,iv),1),AngleTran);
                                            
                        DirMap(abs(AngleTranM-ThetaCMat)<0.5*ThetaStep)=1;  %% needs the be corrected fo matrices always same size
                        
                        DirMap(abs(AngleTranM-ThetaCMat-2*pi)<0.5*ThetaStep)=1; %ADDED SINCE NEEDS TO BE MOD 2PI!! %% needs the be corrected fo matrices always same size

                     end
                    
                     Bp(isnan(Bp))=0;
                     Bp=sparse(repelem(Bp,Ti(i,iv),1));                    
                     RTE=sparse(repmat(RT,nSideEls(i,iv),nSideEls(j,jv)));

                     DirMapE=sparse(repmat(DirMap,nSideEls(i,iv),nSideEls(j,jv)));
                     Bp=RTE.*Bp.*DirMapE;
                    %   

                    TrOpLoc(iBlockA:iBlockB,jBlockA:jBlockB)=StepProd*sparse(Bp);
                   
               end
                

            end
        end
        
    end
    
    for pts=1:PlotSize(j)
       
         FDPosVR(1:nnz(FDPosV(:,pts)),pts,j)=(FDPosV(FDPosV(:,pts)>0,pts));
         LGStoreVR(1:nnz(FDPosV(:,pts)),pts,j)=(LGStoreV(FDPosV(:,pts)>0,pts));  
         FDScaleVR(1:nnz(FDPosV(:,pts)),pts,j)=(FDScaleV(FDPosV(:,pts)>0,pts));%
 
    end
    j
    
end
          

  clear AA ACheck a AngleGam AngleRefIn AngleRefOut AngleRefM AngleTran AngleTranM aii ajj 
  clear BB bii bjj Bp Bpl BpU CL 
  clear DaAj DaBj DaCj DbAj DbBj DbCj DirEdge DirMap DirMapE DirSplit DirVec dLds DSi DSj DSiMat DSjMat DSiMatK DSjMatK dummy
  clear EdgeOK FreeEdges FDPos FDPosV FDScale FDScaleV FDVal FDValV GlobMat GlobMatK GlobMatV iVertp iVertpi Ix Iy Iz KCol KeyFD KRow KVal
  clear La Lb LEdge LGam LGStore LGStoreV Lx Ly Lz MaxEdgei MaxEdgej MinEdgei nSideElsj nSideElsi 
  clear nEdgeHiti nedgeHitj order PListP pv RotX RotY RotZ ScInt sElt sEltList sfactx sfacty sfactz sGam
  clear SR ST StepMat TDNorm ThetaAB ThetaC ThetaCMat ThetaGam ThetaI ThetaIi ThetaLoc ThetaLoci ThetaLocL ThetaLocLi ThetaLocU ThetaLocUi tlist1 TurnAxis TurnAxisU TurnX
  clear UVL UVLp values Xa xai XAi xaiMat xaiV xbj XBj xbjMat xbjV xbiV xcjMat xGam XRA XRAi xRay xRayi XRB XRBi xSplit xSpliti xv xvi XYZ
  clear Ya yai YAi yaiMat yaiV ybj YBj ybjMat ybjV ybiV ycjMat yGam YRA YRAi yRay yRayi YRB YRBi ySplit ySpliti yv yvi 
  clear Za zai ZAi zaiMat zaiV zbj ZBj zbjMat zbjV zbiV zcjMat zGam ZRA ZRAi zRay zRayi ZRB ZRBi zSplit zSpliti zv zvi 

%    %%
stage=3
NM=length(TrOpLoc);
LHS=speye(NM)-TrOpLoc;

VecInf=LHS\ScVec;

clear LHS

stage=4
%%
%Compute final density
FD=zeros(max(PlotSize),nOmega);
IntVal=zeros(max(PlotSize),nOmega);
for j=1:nOmega    
    
        XN=xnodes(Omega(j,1:nVert(j)))';
        YN=ynodes(Omega(j,1:nVert(j)))';
        pv=[XN YN; XN(1) YN(1)]; 
% % % %    
        [points,tri]=distmesh2d(@dpoly,@huniform,PpOm,[min(XN),min(YN); max(XN),max(YN)],pv,pv);
% % % %       
         XP=(points(tri(:,1),1)+points(tri(:,2),1)+points(tri(:,3),1))/3;
         YP=(points(tri(:,1),2)+points(tri(:,2),2)+points(tri(:,3),2))/3;
         ZP=zeros(size(XP));
         
    for pts=1:PlotSize(j)
        
         %Helmholtz version 
        FD(pts,j)=dAbsP(j)*sum(FDScaleVR(1:nnz(FDPosVR(:,pts,j)),pts,j).*VecInf(FDPosVR(1:nnz(FDPosVR(:,pts,j)),pts,j)).*exp(-Adamp(j)*LGStoreVR(1:nnz(FDPosVR(:,pts,j)),pts,j)));
       % Test(pts,j)=sum(VecInf(FDPosVR(1:nnz(FDPosVR(:,pts,j)),pts,j)).*exp(-Adamp(j)*LGStoreVR(1:nnz(FDPosVR(:,pts,j)),pts,j)));
       
       
       if Source==1 && any(j==ScDom)   
             %params.omega*dAbsP(j)
          %Helmholtz
            Zhankre=params.omega*dAbsP(j)*sqrt((Ox-XP(pts))*(Ox-XP(pts))+(Oy-YP(pts))*(Oy-YP(pts))+(Oz-ZP(pts))*(Oz-ZP(pts)));
 		    Zhankim=0.5*Adamp(j)*sqrt((Ox-XP(pts))*(Ox-XP(pts))+(Oy-YP(pts))*(Oy-YP(pts))+(Oz-ZP(pts))*(Oz-ZP(pts)));
% %     %Helmholtz   (params.rhom*params.omega*dAbsP(j)*params.omega*dAbsP(j))*
            Hank=params.rhom*params.omega*params.omega*(abs(besselh(0,1,Zhankre+1i*Zhankim))^2)/16;
          %  Hank=1/(800*pi*pi*sqrt((Ox-XP(pts))*(Ox-XP(pts))+(Oy-YP(pts))*(Oy-YP(pts))+(Oz-ZP(pts))*(Oz-ZP(pts))));
            
           %%integrate H separately over each triangle

%             HankTerm = TriIntegral(@(x,y) (abs(besselh(0,1,params.omega*dAbsP(j)*sqrt((Ox-x).*(Ox-x)+(Oy-y).*(Oy-y))+1i*0.5*Adamp(j)*sqrt((Ox-x).*(Ox-x)+(Oy-y).*(Oy-y)))).^2)/16,points(tri(pts,:),1),points(tri(pts,:),2));       
%              DEASol1(pts)=DEASol1(pts)+Hank;
%              ExactSol1(pts)=ExactSol1(pts)+Hank;

             FD(pts,j)=FD(pts,j)+Hank;
             % FD0(pts,j)=Hank;
          %   IntVal(pts,j)=IntVal(pts,j)+HankTerm;
         end
         
         
         %% Exact solution for rectangle

       ExRect(pts,j)=ExactRect(t0, loss, params.omega, 1./dAbsP(j), 1*sec(t0), params.rhom, XP(pts));

    end
   

%   trisurf(tri,points(:,1),points(:,2),zeros(size(points(:,1))),(FD0(1:PlotSize(j),j)))    
%            hold on
%           axis equal
%           view(2)
%           axis off
if Source==0
  trisurf(tri,points(:,1),points(:,2),zeros(size(points(:,1))),((FD(1:PlotSize(j),j))),'linestyle','none')    
else
  trisurf(tri,points(:,1),points(:,2),zeros(size(points(:,1))),(log10(FD(1:PlotSize(j),j))),'linestyle','none')  
end
    
    hold on
          axis equal
          view(2)
          axis off     


         %figure(2)
       %  subplot(1,2,2) 
%          trisurf(tri,points(:,1),points(:,2),zeros(size(points(:,1))),ExRect(1:PlotSize(j),j))
%           hold on
%           axis equal
%           view(2)

       ErrRect(1:PlotSize(j),j)=(abs(FD(1:PlotSize(j),j)-ExRect(1:PlotSize(j),j)));

end

   if Source==1
       hold on
    plot3(Ox,Oy,Oz,'*')
   end
