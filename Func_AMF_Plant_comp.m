function yout = Func_AMF_Plant_comp(X)
global  Q mup mum Alpha Beta d a Nx Ad_x D_p D_m
Nm = length(X);
N_plant = length(Beta);
N_AMF   = length(Alpha);

yout = zeros(Nm,1);
P = X(1:Nx);
m = X(Nx+1:end);
M = reshape(m,Nx,N_AMF);

% alpha = ALPHA';
% beta  = BETA';

%% Plant
Pi = D_p*Ad_x*P + Q*sum(Alpha.*M,2).*P./(d+P) - ...
    P.*Beta.*sum(M.*a./(a+sum(M,2)-M),2) ...
    - mup*P.*P;

%% AMF
Mj = D_m*Ad_x*M + Beta.*P.*M.*a./(a+(sum(M,2)-M)) - ...
    Alpha.*M.*P./(d+P) - mum*M.*M;

%% Output
yout(1:Nx) = Pi;
yout(Nx+1:end) = Mj(:);
end