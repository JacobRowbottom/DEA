function [rho_exact]=ExactRect(t0, eta, omega, c, L, rho_f, X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exact Solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Domain length L 
%% c=1/|p|, which is speed for Helm and sqrt(omega)/speed for biharmonic

% Loss factor 
%eta = 0.01; 
% \mu
%helmholtz
mu = omega*eta/(2*c); % if input loss
%mu = eta; % if input Adamp not loss
Rho0=1;%rho_f;/(1+(eta^2)/16);
%Rho0=1/(2*c);%*omega


%pause
Xa=X*sec(t0);

% Compute exact solution at X \in [0,L] 
rho_exact = Rho0*(exp(-mu*Xa)+exp(-mu*(2*L-Xa)))./(1-exp(-2*mu*L));
