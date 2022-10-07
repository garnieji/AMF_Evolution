%Write preliminary model for mycorrhiza
clear all
close all

Color = get(gca,'colororder');
Marker = ['o','*','d','^'];

global rH qhp qcp  beta1 beta2 mup mum alpha qcm qhm d ap12 ap21 disp_p disp_m

options = odeset('RelTol',1e-4,'AbsTol',1e-6);

qhp = 3;
qcm = 2;
qcp = 1;
qhm = 1;
Q   = qcm*qhp/(qcp*qhm);
mup = 0.3;
mum = 0.3;
alpha = 5; %0.3; %1.5; 5
beta1 = 0.6;% 0.6; %.6 % .4
beta2 = 0.4; % .4 % .4
Beta = [beta1,beta2];
rH    = 0.0;
d     = 1.2;
ap12  = .1;
ap21  = .1;
disp_p = 0;%.001;
disp_m = 0;%.001;

Ap_inv = [ap21,ap12];

%% Stability of 0
pstar_plus =@(alpha,Beta) alpha*(qhp/qcp+qhm/qcm + qhm/qcm*sqrt((Q-1).^2-4*qhp*mup*mum./(qhm*qcp^2*Beta.^2)))...
       ./(2*(Beta+mup*mum./(qcm*qcp*Beta)))-d;
pstar_minus =@(alpha,Beta) alpha.*(qhp/qcp+qhm/qcm - qhm/qcm*sqrt((Q-1).^2-4*qhp*mup*mum./(qhm*qcp^2*Beta.^2)))...
       ./(2*(Beta+mup*mum./(qcm*qcp*Beta)))-d;
Pstar = pstar_plus(alpha,Beta);
Pstar_plus =@(Beta)(Q+1 + sqrt((Q-1).^2-4*qhp*mup*mum./(qhm*qcp^2*Beta.^2)))...
       ./(2*(1+mup*mum./(qcm*qcp*Beta.^2)));
Pstar_minus =@(Beta) (qhp/qcp+qhm/qcm - qhm/qcm*sqrt((Q-1).^2-4*qhp*mup*mum./(qhm*qcp^2*Beta.^2)))...
       ./(2*(Beta+mup*mum./(qcm*qcp*Beta)));
mstar_plus  =@(alpha,Beta) pstar_plus(alpha,Beta).*(qcm*Beta-qhm*alpha./(d+pstar_plus(alpha,Beta)))/mum; 
mstar_minus =@(alpha,Beta) pstar_minus(alpha,Beta).*(qcm*Beta-qhm*alpha./(d+pstar_minus(alpha,Beta)))/mum; 

Mstar = mstar_plus(alpha,Beta);

Stab0 =  @(alpha,Beta1,Beta2) (Beta2<qhm/qcm*alpha/d.*ap12./(ap12+pstar_plus(alpha,Beta1))) + ( Beta2> qhp/qcp*alpha/d.*ap12./(ap12+pstar_plus(alpha,Beta1)));
stab1 =  @(alpha,Beta1,Beta2) (Beta2-qhm/qcm*alpha/d.*ap12./(ap12+pstar_plus(alpha,Beta1)));
% pstar_minus(alpha*ap12./(ap12+pstar_plus(alpha,Beta1)),Beta2);

dbeta = 1e-2;
betamin = sqrt(4*qhp*mup*mum./(qhm*qcp^2*(Q-1).^2));
betamax = qhp/qcp*alpha/d;

beta = betamin:dbeta:betamax;
[BETA1,BETA2] = meshgrid(beta,beta);
STAB0 = Stab0(alpha,BETA1,BETA2);

figure(12)
clf
hold on
plot(beta,pstar_plus(alpha,beta),'-','color',Color(1,:))
plot(beta,pstar_minus(alpha,beta),'--','color',Color(1,:))
plot(beta,mstar_plus(alpha,beta),'-','color',Color(2,:))
plot(beta,mstar_minus(alpha,beta),'--','color',Color(2,:))

line([betamin,betamax],[0,0],'color','k')
xlim([betamin,betamax])   
% axis([0.4,0.6,0.4,0.6]) 

figure(11)
clf
% plot(beta,Stab0(alpha,11,beta))
contourf(BETA1,BETA2,STAB0,30)
xlabel('$\beta_1$','interpreter','latex','FontSize',16)
ylabel('$\beta_2$','interpreter','latex','FontSize',16)


Dab = @(M) alpha*Ap_inv./(Ap_inv+M) - qcp*d*Beta/qhp;
D = @(M)  qcm*Beta*d./(qhm*alpha*Ap_inv./(Ap_inv+M));

%% Computtation of time evolving solution
Itmax = 5e2;
Stab12 = zeros(1,Itmax);
Stab1  = zeros(1,Itmax);
Stab2  = zeros(1,Itmax);
tstab  = 5e-2;
X0   = zeros(4,Itmax); 
Prop = zeros(2,Itmax);
for it=1:Itmax
   
    %% Initial data
    %%%% Identical
%     p10   = .05*rand;
%     p20   = .05*rand;
%     m0    = .5;%.2*rand+0.4;
%     m10   = m0;%.5; % .2*rand+0.4;
%     m20   = m0;%.5; % .2*rand+0.4;
    
    %%%% Equilibrium
    p10   =  Pstar(1)+.09*rand; 
    p20   = .09*rand; % Pstar(2)+
    m0    = .5;%.2*rand+0.4;
    m10   = Mstar(1)+.2*rand+0.4;% m0; % .5; % .2*rand+0.4;
    m20   = m0;% Mstar(1)+.2*rand+0.4;%.5; % .2*rand+0.4;
    
    X0(:,it) = [p10; p20;m10; m20];
    
    Tfin = 1e2;

    [T,Y] = ode45(@Func_AMF_2Plant_comp_indirect, [0,Tfin], [p10; p20;m10; m20], options);
%     [T1,Y1] = ode45(@Func_AMF_2Plant_comp_indirect, 0:Tfin, [p10; 0;m10; 0], options);
%     [T2,Y2] = ode45(@Func_AMF_2Plant_comp_indirect, 0:Tfin, [0; p20;0; m20], options);
    tstabm1 = mup*Y(end,1)/(Q*alpha/(Y(end,1)+d) -beta1);
    tstabm2 = mup*Y(end,2)/(Q*alpha/(Y(end,2)+d) -beta2);
    prop_P = Y(:,1)./sum(Y(:,1:2),2);
    prop_M = Y(:,3)./sum(Y(:,3:4),2);
    Prop(:,it) = [prop_P(end);prop_M(end)];
    Stab12(1,it) = (1*(prop_P(end)>(1-tstab)) + 2*(prop_P(end)<tstab)).*(Y(end,3)>tstabm1).*(Y(end,4)>tstabm2); % 0= extinction | 1=AMF1 wins | 2= AMF2 wins | 3= AMF1,AMF2 coexist
%     Stab12(1,it) = (1*(Y(end,1)>tstab) + 2*(Y(end,2)>tstab)).*(Y(end,3)>tstabm1).*(Y(end,4)>tstabm2); % 0= extinction | 1=AMF1 wins | 2= AMF2 wins | 3= AMF1,AMF2 coexist
%     Stab12(1,it) = 1*(Y(end,1)>tstab) + 2*(Y(end,2)>tstab); % 0= extinction | 1=AMF1 wins | 2= AMF2 wins | 3= AMF1,AMF2 coexist
%     Stab1(1,it) = 1*(Y1(end,1)>tstab) + 2*(Y1(end,2)>tstab); % 0= extinction | 1=AMF1 wins | 2= AMF2 wins | 3= AMF1,AMF2 coexist
%     Stab2(1,it) = 1*(Y2(end,1)>tstab) + 2*(Y2(end,2)>tstab); % 0= extinction | 1=AMF1 wins | 2= AMF2 wins | 3= AMF1,AMF2 coexist

end


%% Figure

%%% No competition
figure(1)
clf
% for i=0:3
%     is = (Stab12==i);
%     P1is = X0(1,is);
%     P2is = X0(2,is);
%     Mis  = X0(3,is);
%     plot(P1is,P2is,'o','color',Color(i+1,:))
%     hold on
% end
% legend('extinction','P1 wins','P2 wins','coexistence')
scatter(X0(1,:),X0(2,:),25,Prop(1,:),'filled')
xlabel('plant 1')
ylabel('plant 2')
%title('1 plant, X fungi')
set(gca,'fontsize',14)
caxis([0,1])
% 
% %%% Competition between plants
figure(2)
clf
hold on
plot(T,Y(:,1),'g--', T,Y(:,3),'g:','Linewidth',1.5)
plot(T,Y(:,2),'b--', T,Y(:,4),'b:','Linewidth',1.5)
% axis([0 50])
legend({'p1', 'm1', 'p2', 'm2'}, 'Location','best', 'FontSize',12)

% line([0,Tfin],[Pstar(1),Pstar(1)],'color','g')
% line([0,Tfin],[Pstar(2),Pstar(2)],'color','b')
% line([0,Tfin],[Mstar(1),Mstar(1)],'color','g')
% line([0,Tfin],[Mstar(2),Mstar(2)],'color','b')

xlabel('time')
ylabel('Biomass')
%title('1 plant, X fungi')
set(gca,'fontsize',14)




figure(3)
clf
hold on
plot(T,prop_P,'--')
plot(T,prop_M,'-')


