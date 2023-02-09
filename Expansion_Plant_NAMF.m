
%Write preliminary model for mycorrhiza
clear all
close all

options = odeset('RelTol',1e-4,'AbsTol',1e-6);
options_fs = optimoptions('fsolve','FunctionTolerance',1e-12);
options_fmincon = optimoptions('fmincon','FunctionTolerance',1e-12);

%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Time
Tf = 1e2;
dt = 5e-3;
t = 0:dt:Tf; Nt = length(t);

%%% Fonctionnal parameters
global  Q mup mum Alpha Beta d a Nx Ad_x D_p D_m

qhp = 3;
qcm = 2;
qcp = 1;
qhm = 1;
Q   = qcm*qhp/(qcp*qhm);

mup = 0.3; %1/100; %0.3;
mum = 0.3;% 1/20; %0.3;

d     = 1.2;
%% Competition AMF
a = 0.2;
a_ww = a;
a_cw = a;
a_wc = a;
a_cc = a;
%% Equilibrium Multi AMF (P,M) at alpha and beta fixed
Delta_pstar = @(alpha,beta) beta.^2.*(1+Q).^2 ...
                    -4*Q*(1+(length(alpha)-1)./length(alpha).*var(alpha)./mean(alpha).^2).*(beta.^2+mup*mum/length(alpha));
Pstar_plus =@(alpha,beta) (beta.*(1+Q)+sqrt(Delta_pstar(alpha,beta)))./(2*(beta.^2+mup*mum/length(alpha)));
pstar_plus =@(alpha,beta) mean(alpha).*Pstar_plus(alpha,beta)-d;
mstar_plus =@(alpha,beta) pstar_plus(alpha,beta).*(beta-alpha./(mean(alpha)*Pstar_plus(alpha,beta)))/mum;
a_mean_plus     =@(alpha,beta) mean(alpha).*(1-...
            ((length(alpha)-1)./length(alpha).*var(alpha)./mean(alpha).^2)...
            .*1./(beta.*Pstar_plus(alpha,beta)-1));
Pstar_minus =@(alpha,beta) (beta.*(1+Q)-sqrt(Delta_pstar(alpha,beta)))./(2*(beta.^2+mup*mum/length(alpha)));
pstar_minus =@(alpha,beta) mean(alpha).*Pstar_minus(alpha,beta)-d;
mstar_minus =@(alpha,beta) pstar_minus(alpha,beta).*(beta-alpha./(mean(alpha)*Pstar_minus(alpha,beta)))/mum;
a_mean_minus=@(alpha,beta) mean(alpha).*(1-...
            ((length(alpha)-1)./length(alpha).*var(alpha)./mean(alpha).^2)...
            .*1./(beta*Pstar_minus(alpha,beta)-1));       

%% Range Beta = (beta_min,beta_max)
%% AMF trait
alphamin = 0;
alphamax = 1;
N_AMF = 10;
Alpha = linspace(alphamin,alphamax,N_AMF);
% Alpha = alphamax*ones(1,N_AMF);
% Alpha = [alphamin*ones(1,N_AMF/2),alphamax*ones(1,N_AMF/2)];
Nalpha = length(Alpha);
Am_max = alphamax;

Aa = a*(ones(N_AMF,N_AMF)-diag(ones(1,N_AMF)));

%% Plant traits beta
dbeta   = 1e-1;
betamax0 = 50; % 0.6; %.6 % .4
betamax = fsolve(@(Beta) Pstar_plus(Alpha,Beta)*mean(Alpha)-d,betamax0,options_fs);
betamin = sqrt(4*Q*(1+(N_AMF-1)./N_AMF.*var(Alpha)./mean(Alpha).^2)*mup*mum...
       ./(((1+Q).^2-4*Q*(1+(N_AMF-1)./N_AMF.*var(Alpha)./mean(Alpha).^2))*N_AMF))+2*dbeta;

Beta    = betamin + 2*dbeta;
BETA    = betamin:dbeta:betamax;
Nbeta   = length(Beta);

%% Figure theorique A_mean
Amean = a_mean_plus(Alpha,BETA);
figure(4)
clf
yyaxis left
plot(BETA,Amean)
yyaxis right
plot(BETA,pstar_plus(Alpha,BETA))
% plot(BETA,pstar_plus(Alpha,BETA)+Amean)


%% Solving PDE
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
D_m = 0.1;  % diffusion rate
Ap_x  = I_x-dt*D_p*Ad_x;
Am_x  = I_x-dt*D_m*Ad_x;

%% PDE functions
choice = 1;
%% Wiothout competition between AMF
if (choice == 1)
fp = @(p,m,alpha,beta)  Q*sum(alpha.*m,2).*p./(d+p) - p.*beta.*sum(m,2) ...
    - mup*p.*p;
fm = @(p,m,alpha,beta)  beta.*p.*m - ...
    alpha.*m.*p./(d+p) - mum*m.*m;
F = @(X) [Q*sum(Alpha'.*X(2:end)).*X(1)./(d+X(1)) ...
    - X(1).*Beta.*sum(X(2:end)) ...
    - mup.*X(1).*X(1);...
         Beta.*X(1).*X(2:end) - ...
    Alpha'.*X(2:end).*X(1)./(d+X(1)) - mum*X(2:end).*X(2:end)];
elseif (choice == 11)
fp = @(p,m,alpha,beta)  Q*sum(alpha.*m,2).*p./(d+p) - p.*beta.*sum(m,2) ...
    - mup*p.*p;
fm = @(p,m,alpha,beta)  beta.*p.*m - ...
    alpha.*m.*p./(d+p) - mum*m.*sum(m,2);
F = @(X) [Q*sum(Alpha'.*X(2:end)).*X(1)./(d+X(1)) ...
    - X(1).*Beta.*sum(X(2:end)) ...
    - mup.*X(1).*X(1);...
         Beta.*X(1).*X(2:end) - ...
    Alpha'.*X(2:end).*X(1)./(d+X(1)) - mum*X(2:end).*sum(X(2:end))];

%% With competition between AMF
elseif(choice == 2)
fp = @(p,m,alpha,beta)  Q*sum(alpha.*m,2).*p./(d+p) - ...
    p.*beta.*sum(m.*a./(a+sum(m,2)-m),2) ...
    - mup*p.*p;
fm = @(p,m,alpha,beta)  beta.*p.*m.*a./(a+(sum(m,2)-m)) - ...
    alpha.*m.*p./(d+p) - mum*m.*m;

F = @(X) [Q*sum(Alpha'.*X(2:end)).*X(1)./(d+X(1)) ...
    - X(1).*Beta.*sum(X(2:end).*a./(a+sum(X(2:end))-X(2:end))) ...
    - mup.*X(1).*X(1);...
         Beta.*X(1).*X(2:end).*a./(a+(sum(X(2:end))-X(2:end))) - ...
    Alpha'.*X(2:end).*X(1)./(d+X(1)) - mum*X(2:end).*X(2:end)];
end


%% Equilibrium  
Pstar = pstar_plus(Alpha,Beta);
Mstar = mstar_plus(Alpha,Beta);
mstar = max(mstar_minus(Alpha,Beta),0);
Astar_mean = a_mean_plus(Alpha,Beta);
%% Compute equlibrium
        p0 = 20;
        m0 = 10*ones(1,N_AMF);
        Tfe = 1e2;
        [T,X] = ode45(@(t,X) F(X),[0,Tfe],[p0,m0]);

%% Initial conditions
P0 = X(end,1)*(xx<=0)+X(end,1).*exp(-60*(xx)).*(xx>0);
M0 = X(end,2:end)'*(xx<=0)+X(end,2:end)'.*exp(-60*(xx)).*(xx>0);

%% Compute system with ode45
X0 = [P0',M0'];
X0 = X0(:);

if (choice == 1 )
    %% With competition AMF
    [tt,X] = ode45(@(t,y) Func_AMF_Plant_nocomp(y),[0,Tf],X0');
elseif(choice == 2)
    %% With competetion AMF
    [tt,X] = ode45(@(t,y) Func_AMF_Plant_comp(y),[0,Tf],X0');
end
Nt = length(tt);
P = X(:,1:Nx)';
M = reshape(X(:,Nx+1:end),Nt,Nx,N_AMF);
M = permute(M,[2,3,1]);
A_mean = sum(M.*Alpha,2)./(sum(M,2)+(sum(M,2)<=0)).*(sum(M,2)>0);
A_mean = permute(A_mean,[1,3,2]);

%% 
Pbiom = P;
Mbiom = permute(sum(M,2),[1,3,2]);
MM = cumsum(M,2);
PropM = cumsum(M./sum(M,2),2);
Prop_PM = Pbiom./(Mbiom+Pbiom);
%% Figure
figure(12)
clf
for it=  sum(tt<20) % 1:10:Nt-1 %
    clf
    yyaxis left
    hold on
    plot(xx,Pbiom(:,it),'-','LineWidth',3)
%     plot(xx,Mbiom(:,it),'-')
        plot(xx,MM(:,:,it)')
    yyaxis right
    plot(xx,A_mean(:,it))
    plot(xx,PropM(:,:,it))
    plot(xx,Prop_PM(:,it),'-','LineWidth',2)
    line([xx(1),xx(end)],[Astar_mean,Astar_mean])
%     line([xx(1),xx(end)],[A0star_mean,A0star_mean])
% legend(a,[num2str(Alpha(1)),num2str(Alpha(2))])
ylim([0,max(mean(Alpha),1)])
    drawnow
    hold off
end

%%
figure(13)
clf
for it= sum(tt<10)  %1:10:Nt-1
    clf
    yyaxis left
    hold on
    plot(xx,Pbiom(:,it),'-','LineWidth',3)
    plot(xx,Mbiom(:,it),'--')
    ylabel({'Density of plant $P(t,x)$',  'Density of AMF $\displaystyle M(t,x)=\sum_{i=1}^{N_{AMF}}m_i(t,x)$'},'Interpreter', 'latex','FontSize',16)

    yyaxis right
    plot(xx,A_mean(:,it))
    line([xx(1),xx(end)],[Astar_mean,Astar_mean])
ylim([0,max(mean(Alpha),1)])
    drawnow
    hold off
    ylabel({'Mean AMF trait $\displaystyle\sum_{i=1}^{N_{AMF}} \frac{\alpha_i m_i(t,x)}{\sum_{i=1}^{N_{AMF}}m_i(t,x)}$'},'Interpreter', 'latex','FontSize',16)
    xlabel('Space $x$','Interpreter','latex','FontSize',16)
end



%% Figure
figure(2)
clf
for it=1:2:Nt-1
   clf
   yyaxis left
   semilogy(xx,M(:,:,it))
   hold on
   semilogy(xx,P(:,it),'-','LineWidth',3)
      ylim([1e-20,1e2])
   yyaxis right
plot(xx,PropM(:,:,it))

drawnow
end
