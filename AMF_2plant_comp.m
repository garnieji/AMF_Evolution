%Write preliminary model for mycorrhiza
clear all
close all
global rH qhp qcp  beta1 beta2 mup mum alpha qcm qhm d ap12 ap21


qhp = 3;
qcm = 2;
qcp = 1;
qhm = 1;
Q   = qcm*qhp/(qcp*qhm);
mup = 0.3;
mum = 0.3;
alpha = 5; %0.3; %1.5;
beta1 = 0.4; %.6
beta2 = 0.4; % .4
Beta = [beta1,beta2];
rH    = 0.0;
d     = 1.2;
ap12  = .1;
ap21  = .1;


Ap_inv = [ap21,ap12];

pstar =@(alpha,Beta) alpha*(qhp/qcp+qhm/qcm + qhm/qcm*sqrt((Q-1).^2-4*qhp*mup*mum./(qhm*qcp^2*Beta.^2)))...
       ./(2*(Beta+mup*mum./(qcm*qcp*Beta)))-d;
Pstar = pstar(alpha,Beta);
% Pstar2 = alpha*(qhp/qcp+qhm/qcm - qhm/qcm*sqrt((Q-1).^2-4*qhp*mup*mum./(qhm*qcp^2*Beta.^2)))...
%        ./(2*(Beta+mup*mum./(qcm*qcp*Beta)))-d;   
mstar =@(alpha,Beta) pstar(alpha,Beta).*(qcm*Beta-qhm*alpha./(d+pstar(alpha,Beta)))/mum; 
Mstar = mstar(alpha,Beta);

Dab = @(M) alpha*Ap_inv./(Ap_inv+M) - qcp*d*Beta/qhp;
D = @(M)  qcm*Beta*d./(qhm*alpha*Ap_inv./(Ap_inv+M));

%% Initial data
choice_init = 1;
if (choice_init == 1)
    %%%% Identical
    p10   = .02*rand;
    p20   = .02*rand;
    m10   = .5;
    m20   = .5;
elseif (choice_init == 2)
    %%% Close to (P*1,M*1)
    p10   = Pstar(1);  %.2;
    p20   = .002;  %.2;
    m10   = Mstar(1); %.2;
    m20   = .2; %.2;
elseif (choice_init == 3)
    %%% Close to (P*2,M*2)
    p10   = .002;  %.2;
    p20   = Pstar(2);  %.2;
    m10   = .2; %.2;
    m20   = Mstar(2); %.2;
end
%%

Tfin = 1e3;

options = odeset('RelTol',1e-4,'AbsTol',1e-6);
choice_comp = 2 ;
if (choice_comp==1)
%%% Indirect competition between plants through AMF 
[T,Y] = ode45(@Func_AMF_2Plant_comp_indirect, 0:.1:Tfin, [p10; p20;m10; m20], options);
[T1,Y1] = ode45(@Func_AMF_2Plant_comp_indirect, 0:.1:Tfin, [p10; 0;m10; 0], options);
[T2,Y2] = ode45(@Func_AMF_2Plant_comp_indirect, 0:.1:Tfin, [0; p20;0; m20], options);
elseif (choice_comp ==2)
%%% Direct competition between plants
[T,Y] = ode45(@Func_AMF_2Plant_comp_direct, 0:.1:Tfin, [p10; p20;m10; m20], options);
[T1,Y1] = ode45(@Func_AMF_2Plant_comp_direct, 0:.1:Tfin, [p10; 0;m10; 0], options);
[T2,Y2] = ode45(@Func_AMF_2Plant_comp_direct, 0:.1:Tfin, [0; p20;0; m20], options);
end

%%% No competition
figure(1)
clf
hold on
plot(T1,Y1(:,1),'g--', T1,Y1(:,3),'g:','Linewidth',1.5)
plot(T2,Y2(:,2),'b--', T2,Y2(:,4),'b:','Linewidth',1.5)
% axis([0 50])
% legend({'p1', 'p2', 'm1', 'm2'}, 'Location','best', 'FontSize',12)

line([0,Tfin],[Pstar(1),Pstar(1)],'color','g')
line([0,Tfin],[Pstar(2),Pstar(2)],'color','b')
line([0,Tfin],[Mstar(1),Mstar(1)],'color','g')
line([0,Tfin],[Mstar(2),Mstar(2)],'color','b')

xlabel('time')
ylabel('Biomass')
%title('1 plant, X fungi')
set(gca,'fontsize',14)

%%% Competition between plants
figure(2)
clf
hold on
plot(T,Y(:,1),'g--', T,Y(:,3),'g:','Linewidth',1.5)
plot(T,Y(:,2),'b--', T,Y(:,4),'b:','Linewidth',1.5)
% axis([0 50])
legend({'p1', 'm1', 'p2', 'm2'}, 'Location','best', 'FontSize',12)

line([0,Tfin],[Pstar(1),Pstar(1)],'color','g')
line([0,Tfin],[Pstar(2),Pstar(2)],'color','b')
line([0,Tfin],[Mstar(1),Mstar(1)],'color','g')
line([0,Tfin],[Mstar(2),Mstar(2)],'color','b')

xlabel('time')
ylabel('Biomass')
%title('1 plant, X fungi')
set(gca,'fontsize',14)




