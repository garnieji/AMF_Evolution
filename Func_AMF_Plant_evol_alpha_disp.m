function yout = Func_AMF_Plant_evol_alpha_disp(X)
global q_hp q_cm q_hm q_cp beta mup mui d rp ALPHA Ad_alpha dm
Nm = length(X);
yout = zeros(Nm,1);
P = X(1);
M = X(2:end);
alpha = ALPHA';
pnew = (q_hp*rp + sum(alpha.*X(2:end)).*q_hp./(d+X(1)) ...
        -q_cp*sum(beta*X(2:end))...
        -mup*X(1)).*X(1) ;
Mj = dm*Ad_alpha*X(2:end) + X(2:end).*(q_cm*beta*X(1) - ...
      q_hm*alpha.*X(1)./(d+X(1)) - mui*X(2:end));
yout(1) = pnew;
yout(2:end) = Mj;

end