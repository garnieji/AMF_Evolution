clear all
close all
% format LongE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Michorizer model with frequency dependence and evolution of alpha no space   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time parameter
Tf = 300;
%% Parameter of the model
global q_hp q_cm q_hm q_cp mup mui d rp Aa a Ad_alpha ALPHA dm BETA Ap Ad_beta dp dalpha
q_hp  = 3; % q>1
q_cm  = 2;
q_hm  = 1;
q_cp  = 1;
%% Trait beta
beta  = 0.4;
beta_min = 0;
beta_max = 1;
dbeta = 0.1;
% BETA = beta_min:dbeta:beta_max;
BETA = .6;
Nbeta = length(BETA);
% alpha = beta;
%% Trait alpha
alphamin = 0;
alphamax = 10;
dalpha = 0.1;
ALPHA  = alphamin:dalpha:alphamax;
Nalpha = length(ALPHA);

mup = 1/100; % 1/100 %0.3
mui = 1/20; % 1/20 % 0.03

d = 1.2;

rp = 0;
%% Competition terms
a = .2; % 0.1;
ap = 0.1;

N_AMF = Nalpha;
N_plant = Nbeta;
Aa = a*(ones(N_AMF,N_AMF)-diag(ones(1,N_AMF)));
Ap = ap*(ones(N_plant,N_plant)-diag(ones(1,N_plant)));


% Diffusion matrix alpha 
e = ones(Nalpha,1);
I_alpha  = spdiags(e,0,Nalpha,Nalpha);
Ad_alpha = spdiags([e -2*e e],-1:1,Nalpha,Nalpha);
Ad_alpha(1,1) = -1;
Ad_alpha(end,end) = -1;
Ad_alpha = Ad_alpha/(dalpha^2);
% Diffusion matrix beta 
e = ones(Nbeta,1);
I_beta  = spdiags(e,0,Nbeta,Nbeta);
Ad_beta = spdiags([e -2*e e],-1:1,Nbeta,Nbeta);
Ad_beta(1,1) = -1;
Ad_beta(end,end) = -1;
Ad_beta = Ad_beta/(dbeta^2);


dm = 0.1;  % mutation rate AMF
dp = 0;  % mutation rate plant
% %% Space x
% xmin = -10;
% xmax = 10;
% dx = 0.1;
% xx = xmin:dx:xmax;
% Nx = length(xx);
% 
% % Diffusion matrix
% e = ones(Nx,1);
% I  = spdiags(e,0,Nx,Nx);
% Ad_x = spdiags([e -2*e e],-1:1,Nx,Nx);
% Ad_x(1,1) = -1;
% Ad_x(end,end) = -1;
% Ad_x = Ad_x/(dx^2);

% D_p = 0.1;
% D_m = 0.1;  % diffusion rate
%% Functionnal response and interactions MARIA
% sM = @(M) sum(M.*ones(1,N_AMF) - diag(M))'*dalpha;
% Gamma = @(M) Aa*M./(sM(M)+(sM(M)<=0));
% fp = @(alpha,P,M) P.*(q_hp*rp + q_hp*sum(alpha.*M)*dalpha./(d+P) ...
%     -q_cp*beta*sum(M.*Gamma(M)./(Gamma(M) + sM(M) + ((Gamma(M) + sM(M))<=0) ))*dalpha ...
%     -mup*P) ;
% sM = @(M) sum(M.*ones(1,N_AMF) - diag(M))';
% Gamma = @(M) Aa*M./(sM(M)+(sM(M)<=0));
% fp = @(alpha,P,M) P.*(q_hp*rp + q_hp*sum(alpha.*M)./(d+P) ...
%     -q_cp*beta*sum(M.*Gamma(M)./(Gamma(M) + sM(M) + ((Gamma(M) + sM(M))<=0) )) ...
%     -mup*P) ;
% fm = @(alpha,P,M) M.*( (q_cm*beta*Gamma(M)./(Gamma(M) + sM(M) + ((Gamma(M) + sM(M))<=0) )...
%     - q_hm*alpha./(d+P)).*P - mui.*M );

%% Initial data
P0 = 0.1*rand(1,N_plant);
M0 = 0.1*rand(1,N_AMF);
% M0 = 0.1*(ALPHA>1);
t = 0; it = 0; 
tt = t;
dt = 0.1;
Pnew = P0';  PP = P0;  
Mnew = M0'; MM = M0;
%% Explicit sceme
% MM_b = sum(M0);
% MM_d_new = MM./MM_b; MM_d_old = 0;
% while (t<Tf)&&(sum(abs(MM_d_new-MM_d_old))>1e-6)
%     Pold = Pnew; Mold = Mnew; MM_d_old = MM_d_new;
%     Pnew = Pold + dt*fp(ALPHA',Pold,Mold);
%     Mnew = (I -  dt*dm*Ad_alpha)\(Mold + dt*fm(ALPHA',Pold,Mold));
%     MM_d_new = Mnew'./sum(Mnew);
%  
%     it = it + 1; t = t+dt;
%     tt = [tt;t];
%     PP = [PP;Pnew];
%     MM = [MM;Mnew'];
%     MM_bt = sum(Mnew);
%     MM_b = [MM_b;MM_bt];
%     
% end
%% ode45 scheme competition
X0 = [P0,M0];
[t,X] = ode45(@(t,y) Func_AMF_Plant_evol_alpha_comp_continuous(y),[0,Tf],X0);
tt = 0:dt:Tf;
PP = interp1(t,X(:,1),tt);
MM = interp1(t,X(:,2:end),tt);
MM_b = sum(MM,2);

%% ode45 scheme no competition
% X0 = [P0,M0];
% [t,X_nc] = ode45(@(t,y) Func_AMF_Plant_evol_alpha_nodisp(y),[0,Tf],X0);
% tt_nc = 0:dt:Tf;
% PP_nc = interp1(t,X(:,1),tt);
% MM_nc = interp1(t,X(:,2:end),tt);
% MM_b_nc = sum(MM,2);

%% Plot biommass 
% PP_b = sum(PP,2)*dx;

figure(1)
clf
hold on
plot(tt,PP,'--')
plot(tt,MM_b,'-o')
% xlim([0,Tf])
drawnow
hold off

%% Plot of M distribution over time over space trait
MM_d = MM./(MM_b*dalpha);
figure(2)
clf

for i = 1:10:length(tt)
plot(ALPHA,MM_d(end,:))
hold on
plot(ALPHA,MM_d(i,:))
drawnow
hold off
% pause
end
% ylim([0,3e-3])

% for It = 1:10:length(tt)
%     figure(3)
%     clf
%     hold on
%     plot(xx,MM_d(It,:),'-')
%     drawnow
%     pause(0.1)
%     hold off
%     
% %     figure(5)
% %     clf
% %     yyaxis left
% %     Fm = fm(PP(It,:),MM(It,:));
% %     plot(xx,Fm)
% % %     axis([xmin,xmax,0,max(Fm,1)])
% %     yyaxis right
% %     plot(xx,MM(It,:))
% %     axis([xmin,xmax,0,1.01*mc_sstar])
% %     drawnow
% end





