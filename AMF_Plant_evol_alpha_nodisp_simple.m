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
rp = 0;%.2; %.2;
betamin = sqrt(4*Q*mup*mui./((Q-1).^2));
betamax = 13;
%BETA = [linspace(betamin,1,20),linspace(1.1,betamax,20)];
BETA = [logspace(log10(betamin),0,20),logspace(log10(1.1),log10(betamax),20)];
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
Prop_cheater_mutualist = zeros(1,Nbeta);

%% Diffusion matrix
e = ones(Nalpha,1);
I_alpha  = spdiags(e,0,Nalpha,Nalpha);
Ad_alpha = spdiags([e -2*e e],-1:1,Nalpha,Nalpha);
Ad_alpha(1,1) = -1;
Ad_alpha(end,end) = -1;
Ad_alpha = Ad_alpha/(dalpha^2);

dm = 1e-3;  % mutation rate

% ib = 6;
parfor ib = 1:Nbeta
    beta  = max(betamin,BETA(ib));
    if (rp==0)
    Pstar  = (Q+1+sqrt( (Q-1)^2 -4*Q*mup*mui./(beta.^2)) )./(2*(beta+mup*mui./beta));
    alpha_c = max(beta*d/Q,d./Pstar);
    %% Equilibrium (P,M) at alpha and beta fixed    
    C = [1,-d/Pstar,0,-(z0)^3*dm];
    abar_approx = max(abs(roots(C)));
    else
         options = optimoptions('fsolve','MaxFunctionEvaluations',2e3);
        P0star  = (Q+1+sqrt( (Q-1)^2 -4*Q*mup*mui./(beta.^2)) )./(2*(beta+mup*mui./beta));
        C0 = [1,-d/P0star,0,-(z0)^3*dm];
        Abar0 = max(abs(roots(C0)));
        F1pstar = @(X,Y) -Y.*(beta^2+mui*mup).*(X.^3) ...
            + (d.*(beta^2+mui*mup) +beta*(Q+1).*Y - rp*mui).*(X.^2) ...
            -(Y*Q -beta*d*(Q+1)).*X + Q*d;
        F2pstar = @(X,Y) (Y.^2).*(Y.*X-d) - z0^3*dm.*X;
        Fpstar =@(X,Y) [F1pstar(X,Y) ; F2pstar(X,Y)];
        X0 = [P0star+0.4,Abar0];
        X = fsolve(@(X) Fpstar(X(1),X(2)),X0,options);
        Pstar = X(1);
        abar_approx = X(2);
        alpha_c = max(beta*d/Q,d./Pstar);
    end
    
    %% Functionnal response with competition on maintenance cost (strong)
    fp = @(alpha,P,M)  P.*(rp+ (Q*sum(alpha.*M*dalpha) -beta*sum(M*dalpha))-mup*P);
    fm = @(alpha,P,M)  M.*( (beta - alpha).*P - mui.*sum(M).*dalpha);
    F = @(X) [fp(ALPHA',X(1),X(2:end));...
        dm*Ad_alpha*X(2:end)+fm(ALPHA',X(1),X(2:end))];
    %% Trait distribution equilibrium
    p0 = 9;
    m0 = 4*ones(1,N_AMF);
    if (beta<4)
        Tfe = 1e2; % 1e2;
    else
        Tfe = 3e2;
    end
    [T,X] = ode45(@(t,X) F(X),[0,Tfe],[p0,m0]);
    PP(ib) = X(end,1);
    M = X(end,2:end);
    MM(ib,:) = M;
    M_b = sum(M)*dalpha;
    MM_b(ib) = M_b;
    MM_d(ib,:) = M./(M_b);
    Abar(ib) = sum(ALPHA.*M)./sum(M);
%     prop = sum(M.*(ALPHA<=alpha_c))./sum(M);
%     Prop_cheater_mutualist(ib) = prop;
    
    %% Approximations   
    P_approx    = abar_approx*Pstar-d;
    M_b_approx  = P_approx.*(beta-1./Pstar)/mui;
    eta = (P_approx./(P_approx+d)/dm)^(1/3);
    M_d_approx = airy( eta*ALPHA-z0 );
    
    PP_approx(ib) = P_approx;
    MM_b_approx(ib) = M_b_approx;
    MM_d_approx(ib,:) = M_d_approx./sum(M_d_approx*dalpha);
    Abar_approx(ib) = abar_approx;
    
end

%% Figures configurations
% set up figure (Rebecca's version)
% colour scheme
bluegradient2 = [ [97 148 188]/255; [116 162 199]/255; [134 177 210]/255; [152,191,222]/255; ...
    [171 205 233]/255; [190 220 244]/255; [208 234 255]/255];
red2 = [188 70 118]/255;
redgradient = [ red2; [198 98 139]/255; [208 126 160]/255; [218 136 170]/255];

%% Plot of  equilibrium M distribution over space trait
Existence = BETA*d/Q;
Cheat     = BETA*d;
PSTAR  = (Q+1+sqrt( (Q-1)^2 -4*Q*mup*mui./(BETA.^2)) )./(2*(BETA+mup*mui./BETA));
ALPHA_critic = max(Existence,d./PSTAR);
% ALPHA_critic = max(Cheat,d./PSTAR);

[P_max,Ipmax] = max(PP);
Beta_P_max = BETA(Ipmax);
[M_max,Immax] = max(MM_b);
Beta_M_max = BETA(Immax);

figure(1)
clf
% yyaxis left
semilogx(BETA,PP_approx,'--','LineWidth',2,'color', 'k')
hold on
scatter(BETA,PP,'v','filled','MarkerFaceColor','k')
scatter(BETA,MM_b,'filled','MarkerFaceColor', bluegradient2(2,:))
semilogx(BETA,MM_b_approx,'-','LineWidth',2,'Color', bluegradient2(2,:))
ylabel('Plant and total AMF biomass','Interpreter','latex','FontSize',20)
yyaxis right
% semilogx(BETA,Existence,'-','LineWidth',3,'Color',[0.85,0.85,0.85])
semilogx(BETA,ALPHA_critic,'-.','LineWidth',3,'Color',[0.85,0.85,0.85])
hold on
scatter(BETA,Abar,'filled','MarkerFaceColor',red2)
semilogx(BETA,Abar_approx,'-.','LineWidth',2,'Color',red2)

Id = [6,20,31,40];
for im = 1:length(Id)
    id = Id(im);
    line([BETA(id),BETA(id)],[0,Abar(id)],'linestyle',':','color',[0.7,0.7,0.7],'linewidth',1.5)
    plot(BETA(id),Abar(id),'*k', 'MarkerSize', 10,'LineWidth',2)%,'MarkerEdgeColor',redgradient(im,:)); 
end

ax = gca;
ax.YAxis(2).Color = red2;

ylabel({'Mean mutualistic ', 'investment of AMF, $\overline{\alpha}$'},'Interpreter','latex','FontSize',20)
xlim([0,betamax+0.5])
xlabel('Carbon supply rate, $\log(\beta)$','Interpreter','latex','FontSize',20)
% saveas(gca,'fig_biomass_mean_trait_beta_rp0.eps','epsc')


%%
Id = [6,20,31,40];
BB = floor(BETA(Id));
BB(1) = 0.4;
Ia = 1:10:Nalpha;
style = {'-','--','-.',':'};
bluegradient = ["#00008b"; "#1d289b"; "#314fac"; "#5877bd"; "#759fcd"; "#92c6de"; "#afeeee"];
p = zeros(1,length(Id));
dp = cell(1,length(Id));
figure(2)
clf
hold on
for im = 1:length(Id)
    p(im) = scatter(ALPHA(Ia),MM_d(Id(im),Ia),'filled','MarkerFaceColor',redgradient(im,:));%, 'DisplayName', ['\beta=', num2str(BETA(im))]) %,'LineWidth', 1
    dp{im} = ['\beta=', num2str(BB(im))];
    plot(ALPHA,MM_d_approx(Id(im),:),'-','LineWidth',1.5,'color',redgradient(im,:))
end
for im = 1:length(Id)
    ia = sum(ALPHA<Abar(Id(im)));
    line([Abar(Id(im)),Abar(Id(im))],[0,MM_d_approx(Id(im),ia)],'linestyle',':','color',[0.85,0.85,0.85],'linewidth',1.5)
    plot(Abar(Id(im)),MM_d_approx(Id(im),ia),'*k', 'MarkerSize', 10,'LineWidth',2)%,'MarkerEdgeColor',redgradient(im,:)); 
end
legend(p,dp,'FontSize',16)
% line([Abar,Abar],[0.2,1.8])
% line([Abar_approx,Abar_approx],[0.2,1.8],'color','r')
ylim([0,2.5])
ylabel({'Distribution of the mutualistic ','investment in the AMF community'},'Interpreter','latex','FontSize',20)
xlabel('AMF trait, $\alpha$,','Interpreter','latex','FontSize',20)
% saveas(gca,'fig_distrib_alpha_rp0.eps','epsc')

%% Proportion
Prop_cheater = zeros(1,Nbeta);
Prop_approx  = zeros(1,Nbeta);
Proba_cheater_common_ancestor = zeros(1,Nbeta);
Proba_cheater_common_approx   = zeros(1,Nbeta);
Abar_common_ancestor = sum(ALPHA.*(MM.^2),2)./sum(MM.^2,2);
Abar_common_approx   = sum(ALPHA.*(MM_d_approx.^2),2)./sum(MM_d_approx.^2,2);
parfor ib =1:Nbeta
    mm = MM(ib,:);
    ac = real(ALPHA_critic(ib));
    alpha = [linspace(alphamin,ac,300),linspace(ac+dalpha,alphamax,600)];
    mma = interp1(ALPHA,mm,alpha);
    Ib = alpha<=ac;
    Prop_cheater(ib) = trapz(alpha(Ib),mma(Ib))./trapz(alpha,mma);
    Proba_cheater_common_ancestor(ib) = trapz(alpha(Ib),mma(Ib).^2)./trapz(alpha,mma.^2);
    
    mm_approx = MM_d_approx(ib,:)*MM_b_approx(ib);
    mma_approx = interp1(ALPHA,mm_approx,alpha);
    Prop_approx(ib) = trapz(alpha(Ib),mma_approx(Ib))./trapz(alpha,mma_approx);
    Proba_cheater_common_approx(ib) = trapz(alpha(Ib),mma_approx(Ib).^2)./trapz(alpha,mma_approx.^2);    
end
%%
figure(3)
clf
yyaxis left
semilogx(BETA,Prop_approx,'-','color',bluegradient2(1,:))
hold on
scatter(BETA,Prop_cheater,'filled','MarkerFaceColor',bluegradient2(1,:))
ylabel({'Proportion of parasits',' in AMF community'},'Interpreter','latex','FontSize',20)
ylim([0.25,0.7])
yyaxis right
scatter(BETA,PP,'v','filled','MarkerFaceColor', 'k')
semilogx(BETA,PP_approx,'--','LineWidth',2,'color', 'k')

semilogx(BETA,real(ALPHA_critic),'-.','LineWidth',3,'Color',[0.85,0.85,0.85])

scatter(BETA,Abar,'filled','MarkerFaceColor',red2)
semilogx(BETA,Abar_approx,'-.','LineWidth',2,'Color',red2)

ax = gca;
ax.YAxis(2).Color = 'k';
% ax.YAxis(1).Color = bluegradient2(1,:);


ylabel({'Mean trait of AMF community','and plant biomass'},'Interpreter','latex','FontSize',20)
xlim([0,betamax])
xlabel('Carbon supply rate ($\beta$)','Interpreter','latex','FontSize',20)

%% 
Id = 20;%6;  % beta = 0.4 optimal plant biomass
Color = get(gca,'colororder');
Marker = ['o','*','d','^','v'];

mm = MM(Id,:);
ac = real(ALPHA_critic(Id));
alpha = [linspace(alphamin,ac,300),linspace(ac+dalpha,alphamax,600)];
mma = interp1(ALPHA,mm,alpha);
mma = mma./trapz(alpha,mma);

figure(5)
clf
hold on
Ia = alpha<=ac;
aa = [alpha(Ia), fliplr(alpha(Ia))];
inBetween = [mma(Ia), zeros(1,sum(Ia))];
f = fill(aa, inBetween,red2);
set(f,'EdgeColor','none','FaceAlpha', 0.5)

ia = sum(Ia);

Ia = alpha>=ac;
aa = [alpha(Ia), fliplr(alpha(Ia))];
inBetween = [mma(Ia), zeros(1,sum(Ia))];
f = fill(aa, inBetween,bluegradient2(1,:));
set(f,'EdgeColor','none','FaceAlpha', 0.5)

line([ac,ac],[0,mma(ia)],'color','k','linestyle','--','linewidth',2)
plot(alpha,mma,'k','LineWidth',2)
axis([0,2,0,2.5])
ylabel({'Trait distribution in','the AMF community'},'Interpreter','latex','FontSize',20)
xlabel('AMF trait ($\alpha$)','Interpreter','latex','FontSize',20)

%% Fig Common ancestors proba


figure(3)
clf
Color = get(gca,'colororder');
yyaxis left
semilogx(BETA,Proba_cheater_common_approx,'-','color',bluegradient(1,:))
hold on
scatter(BETA,Proba_cheater_common_ancestor,'^','filled','MarkerFaceColor',bluegradient(1,:))

semilogx(BETA,Prop_approx,'-','color',bluegradient2(1,:))
hold on
scatter(BETA,Prop_cheater,'filled','MarkerFaceColor',bluegradient2(1,:))

ylabel({'Probability of parasitic', 'common ancestor \emph{vs}', 'proportion of paraistic AMF'},'Interpreter','latex','FontSize',20)
ylim([0.25,0.9])

yyaxis right

scatter(BETA,Abar,'filled','MarkerFaceColor',red2)
semilogx(BETA,Abar_approx,'-.','LineWidth',2,'Color',red2)

scatter(BETA,Abar_common_ancestor,'d','filled','MarkerFaceColor',Color(2,:))
semilogx(BETA,Abar_common_approx,'-.','LineWidth',2,'Color',Color(2,:))


ax = gca;
ax.YAxis(2).Color = 'k';
ax.YAxis(1).Color = 'k';


ylabel({'Mean trait of AMF community','and common ancestor'},'Interpreter','latex','FontSize',20)
xlim([0,betamax])
xlabel('Carbon supply rate ($\beta$)','Interpreter','latex','FontSize',20)


%% Fig distrib common ancestor
Id = 6;%6;  % beta = 0.4 optimal plant biomass
Marker = ['o','*','d','^','v'];
Ia = 1:10:Nalpha;

mm = MM(Id,:)./sum(MM(Id,:));
mm_approx = MM_d_approx(Id,:)./sum(MM_d_approx(Id,:));
mm_common = MM(Id,:).^2./sum(MM(Id,:).^2); 
mm_common_approx = MM_d_approx(Id,:).^2./sum(MM_d_approx(Id,:).^2);

figure(5)
clf
hold on
scatter(ALPHA(Ia),mm(Ia),'filled','MarkerFacecolor',redgradient(1,:))
plot(ALPHA,mm_approx,'-','color',redgradient(1,:))
scatter(ALPHA(Ia),mm_common(Ia),'^','filled','MarkerFacecolor',bluegradient(1,:))
plot(ALPHA,mm_common_approx,'-','color',bluegradient(1,:))
axis([0,2,0,0.04])
ylabel({'Trait distribution in','the AMF community','common ancestor'},'Interpreter','latex','FontSize',20)
xlabel('AMF trait ($\alpha$)','Interpreter','latex','FontSize',20)



% Ia = ALPHA<=ALPHA_critic(Id);
% aa = [ALPHA(Ia), fliplr(ALPHA(Ia))];
% inBetween = [MM_d(Id,Ia), zeros(1,sum(Ia))];
% f = fill(aa, inBetween,Color(2,:));
% set(f,'EdgeColor','none','FaceAlpha', 0.5)
% 
% ia = sum(Ia);
% 
% Ia = ALPHA>=ALPHA_critic(Id);
% aa = [ALPHA(Ia), fliplr(ALPHA(Ia))];
% inBetween = [MM_d(Id,Ia), zeros(1,sum(Ia))];
% f = fill(aa, inBetween, Color(1,:));
% set(f,'EdgeColor','none','FaceAlpha', 0.5)
% 
% line([ALPHA_critic(Id),ALPHA_critic(Id)],[0,(MM_d(Id,ia+1)+MM_d(Id,ia))./2],'color','k','linestyle','--','linewidth',2)
% plot(ALPHA,MM_d(Id,:),'k','LineWidth',2)
% axis([0,2,0,2.5])


