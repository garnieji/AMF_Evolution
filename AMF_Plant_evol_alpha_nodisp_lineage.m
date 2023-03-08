clear all
close all
% format LongE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Michorizer model with frequency dependence and evolution of alpha no space   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time parameter
Tf = 5e1;
%% Parameter of the model
global beta ALPHA mui d Ad_alpha Nalpha dm dalpha Ni PPstar MMstar
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
alphamax = 2; %5; %max(2*d/Pstar,5);
dalpha = 0.01;
ALPHA  = alphamin:dalpha:alphamax;
Nalpha = length(ALPHA);
N_AMF = Nalpha;


%% Plant
rp = 0; %.02;
betamin = sqrt(4*Q*mup*mui./((Q-1).^2));
betamax = 13;
% BETA = [linspace(betamin,1,20),linspace(1.1,betamax,20)];
BETA = 0.4;
Nbeta = length(BETA);

%% Equilibrium
beta = BETA;
Pstar  = (Q+1+sqrt( (Q-1)^2 -4*Q*mup*mui./(beta.^2)) )./(2*(beta+mup*mui./beta));

%% Diffusion matrix
e = ones(Nalpha,1);
I_alpha  = spdiags(e,0,Nalpha,Nalpha);
Ad_alpha = spdiags([e -2*e e],-1:1,Nalpha,Nalpha);
Ad_alpha(1,1) = -1;
Ad_alpha(end,end) = -1;
Ad_alpha = Ad_alpha/(dalpha^2);

dm = 0.01;  % mutation rate

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
PPstar = X(end,1);
MMstar = X(end,2:end);

%% Lineage
MMi = zeros(Nalpha,Nalpha);
parfor ii =1:Nalpha
    mi0 = zeros(1,Nalpha);
    mi0(ii) = MMstar(ii);
    % mi0 = diag(e); %(ALPHA'<alphamin+dalpha); %
    % mi0 = mi0(:,end);
    % Ni  = size(mi0,2);
    % Mi0 = mi0(:);
    fm_star = @(alpha,M)  M.*( (beta - alpha./(d+PPstar)).*PPstar - mui.*sum(MMstar).*dalpha);
    Fi =@(X) dm*Ad_alpha*X+fm_star(ALPHA',X);
    [ti,Xi] = ode45(@(t,X) Fi(X),[0,Tfe],mi0);
    MMi(ii,:) = Xi(end,:);
end
% MMi = reshape(Xi(end,:),[Nalpha,Ni]);
% I = (ALPHA==1);
% Prop = mean(MMi./MMstar,2)./dalpha;
Prop = MMi./MMstar./dalpha;

%% Plot of  equilibrium M distribution over space trait

figure(1)
clf
plot(ALPHA,cumsum(MMi,1)')
hold on
plot(ALPHA,MMstar,'LineWidth',2)

%% 
Color = get(gca,'colororder');

figure(2)
clf
hold on
plot(ALPHA,MMstar.^2./(sum(MMstar.^2*dalpha)),'-','LineWidth',6,'Color',[0.8,0.8,0.8])
plot(ALPHA,Prop(:,50),'--','Color',Color(1,:),'LineWidth',2)
% scatter(ALPHA(1:3:end),Prop(1:3:end,50),'filled','MarkerFaceColor',Color(1,:))
plot(ALPHA,MMstar,'-','LineWidth',2,'Color',Color(2,:))
ylabel({'AMF distribution','Fixation probability'},'Interpreter','latex','FontSize',20)
xlabel('Carbon supply rate ($\beta$)','Interpreter','latex','FontSize',20)

AAbar = sum(ALPHA.*MMstar.^2)./sum(MMstar.^2);
Abar = sum(ALPHA.*MMstar)./sum(MMstar);


%% Function Fm_star
function [yout] = Fm_star(X)
global beta ALPHA mui d Ad_alpha Nalpha dm dalpha Ni PPstar MMstar
Mi = reshape(X,[Nalpha,Ni]);
fm_star = Mi.*( (beta - ALPHA'./(d+PPstar)).*PPstar - mui.*sum(MMstar).*dalpha);
yi = dm*Ad_alpha*Mi+fm_star;
yout = yi(:);
end


% saveas(gca,'fig_biomass_mean_trait_beta.eps','epsc')
% 
% 
% %%
% Id = [4,20,30,40];
% Ia = 1:10:Nalpha;
% bluegradient = ["#00008b"; "#1d289b"; "#314fac"; "#5877bd"; "#759fcd"; "#92c6de"; "#afeeee"];
% p = zeros(1,length(Id));
% dp = cell(1,length(Id));
% figure(2)
% clf
% hold on
% for im = 1:length(Id)
%     p(im) = scatter(ALPHA(Ia),MM_d(Id(im),Ia),'filled','MarkerFaceColor',bluegradient(im));%, 'DisplayName', ['\beta=', num2str(BETA(im))]) %,'LineWidth', 1
%     dp{im} = ['\beta=', num2str(BETA(Id(im)))];
%     plot(ALPHA,MM_d_approx(Id(im),:),'-','LineWidth',1.5,'color',bluegradient(im))
% end
% legend(p,dp)
% % line([Abar,Abar],[0.2,1.8])
% % line([Abar_approx,Abar_approx],[0.2,1.8],'color','r')
% ylim([0,2.5])
% ylabel('Trait distribution in the AMF community','Interpreter','latex')
% xlabel('AMF trait ($\alpha$)','Interpreter','latex')
% % saveas(gca,'fig_distrib_alpha.eps','epsc')
% 


