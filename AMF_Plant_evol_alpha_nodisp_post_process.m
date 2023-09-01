clear all
close all
% format LongE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Michorizer model with frequency dependence and evolution of alpha no space   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%% Load Data 
DM = {'0001','0005','001'};
Ddm = {'0.001','0.005','0.01'};
Ndm = length(DM);
PPP = zeros(Ndm,Nbeta);
MMM = zeros(Nbeta,Nalpha,Ndm);
MMM_d = zeros(Nbeta,Nalpha,Ndm);
MMM_b = zeros(Ndm,Nbeta);
AAbar = zeros(Ndm,Nbeta);
PPP_approx = zeros(Ndm,Nbeta);
MMM_d_approx = zeros(Nbeta,Nalpha,Ndm);
MMM_b_approx = zeros(Ndm,Nbeta);
AAbar_approx = zeros(Ndm,Nbeta);
PProp_cheater_mutualist = zeros(Ndm,Nbeta);
for idm = 1:Ndm
load(['AMF_Plant_evol_alpha_nodisp_dm',DM{idm},'.mat'])
PPP(idm,:) = PP;
MMM(:,:,idm) = MM;
MMM_d(:,:,idm) = MM_d;
MMM_b(idm,:) = MM_b;
AAbar(idm,:) = Abar;
PPP_approx(idm,:) = PP_approx;
MMM_d_approx(:,:,idm) = MM_d_approx;
MMM_b_approx(idm,:) = MM_b_approx;
AAbar_approx(idm,:) = Abar_approx;
PProp_cheater_mutualist(idm,:) = Prop_cheater_mutualist;
end


%% Figures configurations
% set up figure (Rebecca's version)
% colour scheme
cmp = colormap(gray);
 blackgradient = [cmp(1,:);cmp(end-35,:);cmp(end-20,:)];
% bluegradient2 = [ [97 148 188]/255; [116 162 199]/255; [134 177 210]/255; [152,191,222]/255; ...
%     [171 205 233]/255; [190 220 244]/255; [208 234 255]/255];
bluegradient2 = [ [97 148 188]/255;  [134 177 210]/255; [171 205 233]/255];
red2 = [188 70 118]/255;
redgradient = [ red2; [198 98 139]/255; [208 126 160]/255; [218 136 170]/255];

%% Plot of  equilibrium M distribution over space trait
Existence = BETA*d/Q;
Cheat     = BETA*d;
PSTAR  = (Q+1+sqrt( (Q-1)^2 -4*Q*mup*mui./(BETA.^2)) )./(2*(BETA+mup*mui./BETA));
ALPHA_critic = BETA*d/Q.*((Q+1-sqrt( (Q-1)^2 -4*Q*mup*mui./(BETA.^2)) )./2); %max(Existence,d./PSTAR);
AC = max(Existence,d./PSTAR);
% ALPHA_critic = max(Cheat,d./PSTAR);

[P_max,Ipmax] = max(PPP);
Beta_P_max = BETA(Ipmax);
[M_max,Immax] = max(MMM_b);
Beta_M_max = BETA(Immax);

p = zeros(1,Ndm);
dp = cell(1,Ndm);
figure(1)
clf

yyaxis right
% semilogx(BETA,Existence,'-','LineWidth',3,'Color',[0.85,0.85,0.85])
semilogx(BETA,real(ALPHA_critic),'-','LineWidth',3,'Color',[0.85,0.85,0.85])
hold on

for idm = 1:Ndm
    yyaxis left
semilogx(BETA,PPP_approx(idm,:),'--','LineWidth',2,'color', blackgradient(idm,:))
semilogx(BETA,MMM_b_approx(idm,:),'-','LineWidth',2,'Color', bluegradient2(idm,:))
hold on
scatter(BETA,MMM_b(idm,:),'filled','MarkerFaceColor', bluegradient2(idm,:))

p(idm) = scatter(BETA,PPP(idm,:),'v','filled','MarkerFaceColor', blackgradient(idm,:))
dp{idm} = ['$d_m$=', Ddm{idm}];

ylabel('Plant and total AMF biomass','Interpreter','latex','FontSize',20)
yyaxis right
scatter(BETA,AAbar(idm,:),'filled','MarkerFaceColor',redgradient(idm,:))
% semilogx(BETA,AAbar_approx(idm,:),'-.','LineWidth',2,'Color',red2)
% 
% Id = [6,20,31,40];
% for im = 1:length(Id)
%     id = Id(im);
%     line([BETA(id),BETA(id)],[0,AAbar(id)],'linestyle',':','color',[0.85,0.85,0.85],'linewidth',1.5)
%     plot(BETA(id),AAbar(id),'*k', 'MarkerSize', 10,'LineWidth',2)%,'MarkerEdgeColor',redgradient(im,:)); 
% end
end
hl = legend(p,dp)
set(hl, 'Interpreter','latex', 'FontSize',16)

ax = gca;
ax.YAxis(2).Color = red2;
ax.YAxis(1).Color = [0,0,0];

ylabel({'Mean mutualistic investment', 'of AMF community, $\overline{\alpha}$'},'Interpreter','latex','FontSize',20)
xlim([0,betamax+0.5])
xlabel('Carbon supply rate, $\beta$ (log scale)','Interpreter','latex','FontSize',20)
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
idm = Ndm;
figure(2)
clf
hold on
for im = 1:length(Id)
    p(im) = scatter(ALPHA(Ia),MMM_d(Id(im),Ia,idm),'filled','MarkerFaceColor',redgradient(im,:));%, 'DisplayName', ['\beta=', num2str(BETA(im))]) %,'LineWidth', 1
    dp{im} = ['\beta=', num2str(BB(im))];
    plot(ALPHA,MMM_d_approx(Id(im),:,idm),'-','LineWidth',1.5,'color',redgradient(im,:))
end
for im = 1:length(Id)
    ia = sum(ALPHA<AAbar(idm,Id(im)));
    line([AAbar(idm,Id(im)),AAbar(idm,Id(im))],[0,MMM_d_approx(Id(im),ia,idm)],'linestyle',':','color',[0.85,0.85,0.85],'linewidth',1.5)
    plot(AAbar(idm,Id(im)),MMM_d_approx(Id(im),ia,idm),'*k', 'MarkerSize', 10,'LineWidth',2)%,'MarkerEdgeColor',redgradient(im,:)); 
end
legend(p,dp,'FontSize',16)
% line([Abar,Abar],[0.2,1.8])
% line([Abar_approx,Abar_approx],[0.2,1.8],'color','r')
ylim([0,2.5])
ylabel({'Distribution of mutualistic','investment in the AMF community'},'Interpreter','latex','FontSize',20)
xlabel('AMF trait, $\alpha$','Interpreter','latex','FontSize',20)
% saveas(gca,'fig_distrib_alpha_rp0.eps','epsc')

%% Proportion
Prop_cheater = zeros(Ndm,Nbeta);
Prop_approx  = zeros(Ndm,Nbeta);
Proba_cheater_common_ancestor = zeros(Ndm,Nbeta);
Proba_cheater_common_approx   = zeros(1,Nbeta);
Abar_common_ancestor = sum(ALPHA.*(MMM.^2),2)./sum(MMM.^2,2);
Abar_common_approx   = sum(ALPHA.*(MMM_d_approx.^2),2)./sum(MMM_d_approx.^2,2);
for idm = 1:Ndm
parfor ib =1:Nbeta
    mm = MMM(ib,:,idm);
    ac = real(ALPHA_critic(ib));
    alpha = [linspace(alphamin,ac,300),linspace(ac+dalpha,alphamax,600)];
    mma = interp1(ALPHA,mm,alpha);
    Ib = alpha<=ac;
    Prop_cheater(idm,ib) = trapz(alpha(Ib),mma(Ib))./trapz(alpha,mma);
    Proba_cheater_common_ancestor(idm,ib) = trapz(alpha(Ib),mma(Ib).^2)./trapz(alpha,mma.^2);
    
    mm_approx = MMM_d_approx(ib,:,idm)*MMM_b_approx(idm,ib);
    mma_approx = interp1(ALPHA,mm_approx,alpha);
    Prop_approx(idm,ib) = trapz(alpha(Ib),mma_approx(Ib))./trapz(alpha,mma_approx);
    Proba_cheater_common_approx(idm,ib) = trapz(alpha(Ib),mma_approx(Ib).^2)./trapz(alpha,mma_approx.^2);    
end
end
%%

idm = Ndm;
figure(3)
clf
yyaxis left
 semilogx(BETA,Prop_approx(idm,:),'-','color',bluegradient2(1,:))
hold on
scatter(BETA,Prop_cheater(idm,:),'filled','MarkerFaceColor',bluegradient2(1,:))
ylabel({'Proportion of parasitic ','AMF in the community'},'Interpreter','latex','FontSize',20)
ylim([0.25,0.7])
yyaxis right
scatter(BETA,PPP(idm,:),'v','filled','MarkerFaceColor', 'k')
semilogx(BETA,PPP_approx(idm,:),'--','LineWidth',2,'color', 'k')

semilogx(BETA,real(ALPHA_critic),'-.','LineWidth',3,'Color',[0.7,0.7,0.7])

scatter(BETA,AAbar(idm,:),'filled','MarkerFaceColor',red2)
semilogx(BETA,AAbar_approx(idm,:),'-.','LineWidth',2,'Color',red2)

ax = gca;
ax.YAxis(2).Color = 'k';
% ax.YAxis(1).Color = 'k';


ylabel({'Plant biomass and mean mutualistic', 'investment of AMF community'},'Interpreter','latex','FontSize',20)
xlim([0,betamax])
xlabel('Carbon supply rate, $\beta$ (log scale)','Interpreter','latex','FontSize',20)

%% 
idm = Ndm;
Id = 20;%6;  % beta = 0.4 optimal plant biomass
Color = get(gca,'colororder');
Marker = ['o','*','d','^','v'];

mm = MMM(Id,:,idm);
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
Abar_common_ancestor = permute(Abar_common_ancestor,[3,1,2]);
Abar_common_approx = permute(Abar_common_approx,[3,1,2]);
%%
idm = Ndm;
p = zeros(1,4);
dp = cell(1,4);
figure(3)
clf
Color = get(gca,'colororder');
yyaxis left
semilogx(BETA,Proba_cheater_common_approx(idm,:),'-','color',bluegradient(1,:))
hold on
p(1) = scatter(BETA,Proba_cheater_common_ancestor(idm,:),'^','filled','MarkerFaceColor',bluegradient(1,:))
dp{1} = 'Proba parasitic CA';
semilogx(BETA,Prop_approx(idm,:),'-','color',bluegradient2(1,:))
hold on
p(2) = scatter(BETA,Prop_cheater(idm,:),'filled','MarkerFaceColor',bluegradient2(1,:))
dp{2} = 'Proportion parasitic AMF';
ylabel({'Probability of parasitic common ancestor', '\emph{vs} proportion of paraistic in AMF community'},'Interpreter','latex','FontSize',20)
ylim([0.25,0.9])

yyaxis right
semilogx(BETA,real(ALPHA_critic),'-','LineWidth',4,'Color',[0.7,0.7,0.7])
%hold on
p(3) = scatter(BETA,AAbar(idm,:),'filled','MarkerFaceColor',red2)
dp{3} = '$\overline{\alpha}$ AMF community';
semilogx(BETA,AAbar_approx(idm,:),'-.','LineWidth',2,'Color',red2)
p(4) = scatter(BETA,Abar_common_ancestor(idm,:),'d','filled','MarkerFaceColor',Color(2,:))
dp{4} = '$\overline{\alpha}$ common ancestor (CA)';
% semilogx(BETA,real(ALPHA_critic),'--','LineWidth',2,'Color','k')
semilogx(BETA,Abar_common_approx(idm,:),'-.','LineWidth',2,'Color',Color(2,:))

ax = gca;
ax.YAxis(2).Color = 'k';
ax.YAxis(1).Color = 'k';

hl = legend(p([3,4,1,2]),dp([3,4,1,2]),'location','northwest')
set(hl, 'Interpreter','latex', 'FontSize',16)

ylabel({'Mean mutualistic investment of ', 'AMF community and common ancestor'},'Interpreter','latex','FontSize',20)
xlim([0,betamax])
xlabel('Carbon supply rate, $\beta$ (log scale)','Interpreter','latex','FontSize',20)
%%


idm = Ndm;
figure(30)
clf
Color = get(gca,'colororder');
yyaxis left
semilogx(BETA,Proba_cheater_common_approx(idm,:),'-','color',bluegradient(1,:))
hold on
scatter(BETA,Proba_cheater_common_ancestor(idm,:),30*Prop_cheater(idm,:),'^','filled')
% CT = cbrewer('seq', 'Blues', 8);
% colormap(CT(3:end,:));
% colorbar('westoutside')
ylim([0.25,0.9])

% semilogx(BETA,Prop_approx(idm,:),'-','color',bluegradient2(1,:))
% hold on
% scatter(BETA,Prop_cheater(idm,:),'filled','MarkerFaceColor',bluegradient2(1,:))

ylabel({'Probability of parasitic','common ancestor'},'Interpreter','latex','FontSize',20)
% xlabel('Proportion of parasit in AMF community','Interpreter','latex','FontSize',20)


yyaxis right
semilogx(BETA,real(ALPHA_critic),'-','LineWidth',4,'Color',[0.7,0.7,0.7])
%hold on
% scatter(BETA,AAbar(idm,:),30,,'filled','MarkerFaceColor',red2)
% semilogx(BETA,AAbar_approx(idm,:),'-.','LineWidth',2,'Color',red2)
scatter(BETA,Abar_common_ancestor(idm,:),40*AAbar(idm,:),'d','filled')
semilogx(BETA,Abar_common_approx(idm,:),'-.','LineWidth',2,'Color',Color(2,:))
% CTR = cbrewer('seq', 'Reds', 8);
% colormap(CTR(3:end,:));
% colorbar('eastoutside')

ax = gca;
ax.YAxis(2).Color = 'k';
ax.YAxis(1).Color = 'k';


ylabel({'Mean mutualistic investment of ', 'AMF community and common ancestor'},'Interpreter','latex','FontSize',20)
xlim([0,betamax])
xlabel('Carbon supply rate, $\log(\beta)$','Interpreter','latex','FontSize',20)

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
ylabel({'Distribution of mutulaistic investment','in AMF community and common ancestor'},'Interpreter','latex','FontSize',20)
xlabel('AMF trait, $\alpha$','Interpreter','latex','FontSize',20)



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


