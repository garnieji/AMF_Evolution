clear all
close all
% format LongE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Michorizer model with frequency dependence and evolution of alpha no space   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time parameter
TF = [70,150];
%% Parameter of the model
% global q_hp q_cm q_hm q_cp beta mup mui d rp
q_hp  = 3; % q>1
q_cm  = 2;
q_hm  = 1;
q_cp  = 1;
Q = q_hp*q_cm/(q_cp*q_hm);

mup = 0.3; % 1/100
mui = 0.3; % 1/20

d = 1.2;                            % plant saturation parameter
%% Plant
rp = 0; %.02;
betamin = sqrt(4*Q*mup*mui./((Q-1).^2));
betamax = 13;
BETA = [linspace(betamin,2,25),linspace(2.1,betamax,15)];
Nbeta = length(BETA);
% alpha = beta;
%% Trait alpha
alphamin = 0;
alphamax = 5;
dalpha = 0.01;                      % stepsize in alpha
ALPHA  = alphamin:dalpha:alphamax;  % vector of alpha values
Nalpha = length(ALPHA);
N_AMF = Nalpha;

% Diffusion matrix
e = ones(Nalpha,1);
I_alpha  = spdiags(e,0,Nalpha,Nalpha);
Ad_alpha = spdiags([e -2*e e],-1:1,Nalpha,Nalpha);
Ad_alpha(1,1) = -1;
Ad_alpha(end,end) = -1;
Ad_alpha = Ad_alpha/(dalpha^2);

dm = 0.01;  % mutation rate

%% Space x
xmin = -5;
xmax = 20;
dx = 0.05;
xx = xmin:dx:xmax;
Nx = length(xx);

% Diffusion matrix
e = ones(Nx,1);
I_x  = spdiags(e,0,Nx,Nx);
Ad_x = spdiags([e -2*e e],-1:1,Nx,Nx);
Ad_x(1,1) = -1;
Ad_x(end,end) = -1;
Ad_x = Ad_x/(dx^2);

D_p = 0.1;
D_m = 0.01;  % diffusion rate

%% Spreading speed / Equilibrium
C_simu = zeros(2,Nbeta);
PSTAR_approx = zeros(1,Nbeta); MSTAR_approx = zeros(1,Nbeta);ABARSTAR_approx = zeros(1,Nbeta);
PSTAR = zeros(1,Nbeta); MSTAR = zeros(1,Nbeta);ABARSTAR = zeros(1,Nbeta);
parfor ib = 1:Nbeta
    beta  = max(betamin,BETA(ib));
%     beta = 8;
if (beta <2)
    Tf = TF(1);
else
    Tf = TF(2);
end
    %% Functional response and interactions MARIA
    fp = @(alpha,P,M)  P.*( q_hp*rp + (Q*sum(ALPHA'.*M*dalpha)./(d+P) -beta*sum(M*dalpha))-mup*P);
    fm = @(alpha,P,M)  M.*( (beta - alpha./(d+P)).*P - mui.*sum(M).*dalpha);
    
    %% Equilibrium
    Pstar  = (Q+1+sqrt( (Q-1)^2 -4*Q*mup*mui./(beta.^2)) )./(2*(beta+mup*mui./beta));
    z0  = -fzero(@(x) airy(1,x),0);
    PA = [1,-d/Pstar,0,-(z0)^3*dm];
    abar_approx = max(abs(roots(PA)));
    pstar    = abar_approx*Pstar-d;
    mstar = pstar.*(beta-1./Pstar)/mui;
    eta = (pstar./(pstar+d)/dm)^(1/3);
    m_d_star = airy( eta*ALPHA'-z0 );
    Mstar = m_d_star./sum(m_d_star*dalpha)*mstar;  
    PSTAR_approx(ib) = pstar; MSTAR_approx(ib) = mstar; ABARSTAR_approx(ib) = abar_approx;
    %% Initial data
    P0 = pstar*(xx<=0);
    M0 = Mstar.*(xx<=0);
    
    t = 0; it = 1;
    dt = 0.01;
    tt = 0:dt:Tf;
    Nt = length(tt);
    Pnew = P0;  Mnew = M0;
    
    Xt = zeros(2,Nt);
    while (it<Nt)
        Pold = Pnew;Mold = Mnew;
        
        Pp = Pold + dt*fp(ALPHA',Pold,Mold);
        Bp = I_x-dt*D_p*Ad_x;
        Ppnew = Bp\(Pp)';
        Pnew = Ppnew';
        
        Mm = Mold + dt*fm(ALPHA',Pold,Mold);
        Am = I_alpha-dt*dm*Ad_alpha;
        Mnew = Am\(dt*D_m*Mold*Ad_x+Mm);
        
       
        it = it + 1; t = t+dt;
        xpt = max(sum(Pnew>(pstar/2)),1);
        xmt = max(sum(sum(Mnew*dalpha,1)>(mstar/2)),1);
        Xt(1,it) = xx(xpt);
        Xt(2,it) = xx(xmt)
    end
        Abar = ALPHA*Mnew./sum(Mnew,1);
        PSTAR(ib) = mean(Pnew(1:100));
        MSTAR(ib) = mean(sum(Mnew(:,1:100),1)*dalpha); 
        ABARSTAR(ib) = mean(Abar(1:100));
%     pxt = polyfit(tt(end-500:end),Xt(end-500:end),1);
%     c_simu = pxt(1);
    c_simu = mean(Xt(:,end-500:end)/(tt(end-500:end)),2);
    C_simu(:,ib) = c_simu;
end

%% Figure speed
% for ib = 1:Nbeta
%     beta  = max(betamin,BETA(ib));
% %     beta = 8;
% if (beta <2)
%     Tf = TF(1);
% else
%     Tf = TF(2);
% end
%     %% Functional response and interactions MARIA
%     fp = @(alpha,P,M)  P.*( q_hp*rp + (Q*sum(ALPHA'.*M*dalpha)./(d+P) -beta*sum(M*dalpha))-mup*P);
%     fm = @(alpha,P,M)  M.*( (beta - alpha./(d+P)).*P - mui.*sum(M).*dalpha);
%     
%     %% Equilibrium
%     Pstar  = (Q+1+sqrt( (Q-1)^2 -4*Q*mup*mui./(beta.^2)) )./(2*(beta+mup*mui./beta));
%     z0  = -fzero(@(x) airy(1,x),0);
%     PA = [1,-d/Pstar,0,-(z0)^3*dm];
%     abar_approx = max(abs(roots(PA)));
%     pstar    = abar_approx*Pstar-d;
%     mstar = pstar.*(beta-1./Pstar)/mui;
%     eta = (pstar./(pstar+d)/dm)^(1/3);
%     m_d_star = airy( eta*ALPHA'-z0 );
%     Mstar = m_d_star./sum(m_d_star*dalpha)*mstar;
%     
%     PSTAR(ib) = pstar; MSTAR(ib) = mstar; ABARSTAR(ib) = abar_approx;
% end
%%
alpha_mean = (alphamin+alphamax)/2;
% df20 = - mui);
%% 1st approx
rd_minus = sqrt((Q*ABARSTAR/d-BETA).*(BETA-ABARSTAR/d)) -mup;
c1_star = sqrt(PSTAR.*D_p.*rd_minus/2);
%% 2nd approx
u_star = (BETA-ABARSTAR/d + mup)./(Q*ABARSTAR/d-BETA+mui);
rd = (Q*ABARSTAR/d-BETA).*u_star-mup;
c2_star = sqrt(PSTAR.*D_p.*rd/2);
%% 3rd
% rd_plus = sqrt((Q*alpha_mean/d-BETA).*(BETA-alpha_mean/d)) -mup;
% c3_star = sqrt(D_p*rd_plus/2);
%% 4th approx
% u_star = (BETA-alpha_mean/d + mup)./(Q*alpha_mean/d-BETA+mui);
% rd = (Q*alpha_mean/d-BETA).*u_star-mup;
% c4_star = sqrt(D_p*rd/2);
%% 5 th approx
% G = [Q*ABARSTAR/d-BETA; mui*MSTAR./PSTAR+mup;-MSTAR./PSTAR.*(BETA-ABARSTAR/d)]';
% GAM = zeros(1,Nbeta);
% for ib = 1:Nbeta
%     GAM(ib) = max(roots(G(ib,:)));
% end    
% rd_5 = (Q*ABARSTAR/d-BETA).*GAM-mup;
G = [Q*ABARSTAR_approx/d-BETA;...
     mui*D_m/D_p*MSTAR_approx./PSTAR_approx+mup;...
     -D_m/D_p*MSTAR_approx./PSTAR_approx.*(BETA-ABARSTAR_approx/d)]';
GAM = zeros(1,Nbeta);
for ib = 1:Nbeta
    GAM(ib) = max(roots(G(ib,:)));
end    
rd_5 = (Q*ABARSTAR_approx/d-BETA).*GAM-mup;
c5_star = sqrt(D_p.*PSTAR_approx.*rd_5/2);

G = [Q*alpha_mean/d-BETA;...
     mui*MSTAR_approx./PSTAR_approx+mup;...
     -MSTAR_approx./PSTAR_approx.*(BETA-alpha_mean/d)]';
GAM = zeros(1,Nbeta);
for ib = 1:Nbeta
    GAM(ib) = max(roots(G(ib,:)));
end    
rd_6 = (Q*alpha_mean/d-BETA).*GAM-mup;
c6_star = sqrt(D_p.*PSTAR_approx.*rd_6/2);

Color = get(gca,'colororder');

figure(11)
clf
% yyaxis left
scatter(BETA,C_simu(1,:),'filled')
hold on
% plot(BETA,c1_star,'-.')
% plot(BETA,c2_star,'--')
% plot(BETA,c3_star,':')
% plot(BETA,c4_star,'d')
plot(BETA,c5_star,'--','Color',Color(1,:))
% plot(BETA,c6_star,'*')
ylabel('Spreading speed (AMF and plant)','Interpreter','latex','FontSize',20)
hold on
% scatter(BETA,C_simu(2,:),'filled')
% line([BETA(1),BETA(end)],[0,0])
ylim([0,0.15])
% yyaxis right
% hold on
% plot(BETA,sqrt(PSTAR))
% % plot(BETA,MSTAR)
% ylabel('Plant and AMF biomass','Interpreter','latex')
xlabel('Carbon suupply rate ($\beta$)','Interpreter','latex','FontSize',20)

figure(2)
clf
yyaxis left
scatter(BETA,ABARSTAR,'filled')
hold on
plot(BETA,ABARSTAR_approx)
yyaxis right
hold on
scatter(BETA,PSTAR,'d','filled')
plot(BETA,PSTAR_approx)
scatter(BETA,MSTAR,'filled')
plot(BETA,MSTAR_approx)




