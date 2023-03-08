clear all
close all
% format LongE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Michorizer model with frequency dependence and evolution of alpha no space   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time parameter
Tf = 5e1;
%% Parameter of the model
% global q_hp q_cm q_hm q_cp beta ALPHA mup mui d rp Ad_alpha Nalpha dm dalpha
q_hp  = 3; % q>1
q_cm  = 2;
q_hm  = 1;
q_cp  = 1;
Q   = q_cm*q_hp/(q_cp*q_hm);

mup = 0.3; % 1/100
mui = 0.3; % 1/20

d = 1.2;

%% Trait alpha
alphamin = 0;
alphamax = 5; %max(2*d/Pstar,5);
dalpha = 0.01;
ALPHA  = alphamin:dalpha:alphamax;
Nalpha = length(ALPHA);
N_AMF = Nalpha;


%% Plant
rp = 0; %.02;
betamin = sqrt(4*Q*mup*mui./((Q-1).^2));
betamax = 13;
BETA = [linspace(betamin,1,20),linspace(1.1,betamax,20)];
Nbeta = length(BETA);

%% Equilibrium
z0  = -fzero(@(x) airy(1,x),0);

PP = zeros(1,Nbeta);
MM = zeros(Nbeta,Nalpha);
MM_d = zeros(Nbeta,Nalpha);
MM_b = zeros(1,Nbeta);
Abar = zeros(1,Nbeta);
PP_approx = zeros(1,Nbeta);
MM_d_approx = zeros(Nbeta,Nalpha);
MM_b_approx = zeros(1,Nbeta);
Abar_approx = zeros(1,Nbeta);

%% Diffusion matrix
e = ones(Nalpha,1);
I_alpha  = spdiags(e,0,Nalpha,Nalpha);
Ad_alpha = spdiags([e -2*e e],-1:1,Nalpha,Nalpha);
Ad_alpha(1,1) = -1;
Ad_alpha(end,end) = -1;
Ad_alpha = Ad_alpha/(dalpha^2);

dm = 0.01;  % mutation rate

% ib = 30;
parfor ib = 1:Nbeta
    beta  = max(betamin,BETA(ib));
    Pstar  = (Q+1+sqrt( (Q-1)^2 -4*Q*mup*mui./(beta.^2)) )./(2*(beta+mup*mui./beta));
    
    
    
    %% Functionnal response with competition on maintenance cost (strong)
    fp = @(alpha,P,M)  P.*( (Q*sum(alpha.*M*dalpha)./(d+P) -beta*sum(M*dalpha))-mup*P);
    fm = @(alpha,P,M)  M.*( (beta - alpha./(d+P)).*P - mui.*sum(M).*dalpha);
    F = @(X) [fp(ALPHA',X(1),X(2:end));...
        dm*Ad_alpha*X(2:end)+fm(ALPHA',X(1),X(2:end))];
    %% Trait distribution equilibrium
    p0 = 9;
    m0 = 4*ones(1,N_AMF);
    
    Tfe = 1e2;
    [T,X] = ode45(@(t,X) F(X),[0,Tfe],[p0,m0]);
    PP(ib) = X(end,1);
    M = X(end,2:end);
    MM(ib,:) = M;
    M_b = sum(M)*dalpha;
    MM_b(ib) = M_b;
    MM_d(ib,:) = M./(M_b);
    Abar(ib) = sum(ALPHA.*M)./sum(M);
    %% Approximations
    %% Equilibrium (P,M) at alpha and beta fixed
    
    C = [1,-d/Pstar,0,-(z0)^3*dm];
    
    abar_approx = max(abs(roots(C)));
    P_approx    = abar_approx*Pstar-d;
    M_b_approx  = P_approx.*(beta-1./Pstar)/mui;
    eta = (P_approx./(P_approx+d)/dm)^(1/3);
    M_d_approx = airy( eta*ALPHA-z0 );
    
    PP_approx(ib) = P_approx;
    MM_b_approx(ib) = M_b_approx;
    MM_d_approx(ib,:) = M_d_approx./sum(M_d_approx*dalpha);
    Abar_approx(ib) = abar_approx;
    
end

%% Plot of  equilibrium M distribution over space trait
Existence = BETA*d/Q;

figure(1)
clf
yyaxis left
scatter(BETA,PP,'d','filled')
hold on
scatter(BETA,MM_b,'filled')
plot(BETA,PP_approx,'--','LineWidth',1)
plot(BETA,MM_b_approx,'-','LineWidth',1)
ylabel('Plant and AMF biomasses','Interpreter','latex')
yyaxis right
hold on
plot(BETA,Existence,'k-','LineWidth',2)
scatter(BETA,Abar,'filled')
plot(BETA,Abar_approx,'-','LineWidth',1)

ylabel('Mean trait $\overline{\alpha}$ of AMF community','Interpreter','latex')
xlim([0,betamax])
xlabel('Carbon supply rate ($\beta$)','Interpreter','latex')
% saveas(gca,'fig_biomass_mean_trait_beta.eps','epsc')


%%
Id = [4,20,30,40];
Ia = 1:10:Nalpha;
bluegradient = ["#00008b"; "#1d289b"; "#314fac"; "#5877bd"; "#759fcd"; "#92c6de"; "#afeeee"];
p = zeros(1,length(Id));
dp = cell(1,length(Id));
figure(2)
clf
hold on
for im = 1:length(Id)
    p(im) = scatter(ALPHA(Ia),MM_d(Id(im),Ia),'filled','MarkerFaceColor',bluegradient(im));%, 'DisplayName', ['\beta=', num2str(BETA(im))]) %,'LineWidth', 1
    dp{im} = ['\beta=', num2str(BETA(Id(im)))];
    plot(ALPHA,MM_d_approx(Id(im),:),'-','LineWidth',1.5,'color',bluegradient(im))
end
legend(p,dp)
% line([Abar,Abar],[0.2,1.8])
% line([Abar_approx,Abar_approx],[0.2,1.8],'color','r')
ylim([0,2.5])
ylabel('Trait distribution in the AMF community','Interpreter','latex')
xlabel('AMF trait ($\alpha$)','Interpreter','latex')
% saveas(gca,'fig_distrib_alpha.eps','epsc')



