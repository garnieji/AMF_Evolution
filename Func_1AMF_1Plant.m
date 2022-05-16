function yout = Func_1AMF_1Plant(X)
global  qhp qcp mup mum  qcm qhm d ALPHA BETA
yout = zeros(2,1);
P = X(1);
M = X(2);
    
% alpha = ALPHA;
% beta  = BETA;

%% AMF + Plant
Pi = (ALPHA.*M.*qhp./(d+P) ...
        -qcp*BETA*M...
        -mup*P).*P ;
    
Mj =  M.*(qcm*BETA*P - ...
      qhm*ALPHA.*P./(d+P) - mum*M);
yout(1) = Pi;
yout(2) = Mj;
end