function yout = Func_AMF_Plant_evol_alpha_disp_comp_strong(X)
global q_hp q_cm q_hm q_cp beta ALPHA mup mui d rp a Ad_x Ad_alpha Nx D_p D_m dm Nalpha dalpha
Nm = length(X);
yout = zeros(Nm,1);
P = X(1:Nx);
M = reshape(X(Nx+1:end),[Nx,Nalpha]);

alpha = ALPHA;

Pnew = D_p*Ad_x*P +  (q_hp*rp + sum(alpha.*M.*a./(a+sum(M*dalpha,2))*dalpha,2).*q_hp./(d+P) ...
        -q_cp*sum(beta*M*dalpha,2)...
        -mup*P).*P ;
Mj = D_m*Ad_x*M + dm*M*Ad_alpha + M.*(q_cm*beta*P.*a./(a+sum(M*dalpha,2)) - ...
      q_hm*alpha.*P./(d+P) - mui*sum(M*dalpha,2));
yout(1:Nx) = Pnew;
yout(Nx+1:end) = Mj(:);

end