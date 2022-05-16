%Write preliminary model for mycorrhiza
clear all
close all
global rH qhp qcp mup mum  qcm qhm d ALPHA BETA
qhp = 3;
qcm = 2;
qcp = 1;
qhm = 1;
Q   = qcm*qhp/(qcp*qhm);
mup = 0.3;
mum = 0.3;
rH    = 0.0;
d     = 1.2;

alpha = 2; %0.3; %1.5;
ALPHA = alpha;
beta1 = 0.6;
beta2 = 0.4;
Beta = [beta1,beta2];
BETA = beta2;


%% Growth functions
Fp = @(P,M) P.*(ALPHA.*M.*qhp./(d+P) - qcp*BETA*M  -mup*P) ;
Fm = @(P,M) M.*(qcm*BETA*P - qhm*ALPHA.*P./(d+P) - mum*M);


%% Equilibrium
%%% Jacobian
J = @(P,M) [-P.*ALPHA.*M.*qhp./((d+P).^2)-mup.*P , ...
            P.*(ALPHA.*qhp./(d+P)-qcp*BETA) ; ...
            M.*(qcm.*BETA - qhm.*ALPHA./(d+P)+qhm.*ALPHA.*P./((d+P).^2)),...
            -mum.*M];
%%% Stable
pstar =@(alpha,Beta) alpha*(qhp/qcp+qhm/qcm + qhm/qcm*sqrt((Q-1).^2-4*qhp*mup*mum./(qhm*qcp^2*Beta.^2)))...
       ./(2*(Beta+mup*mum./(qcm*qcp*Beta)))-d;
Pstar = pstar(ALPHA,BETA);
% Pstar2 = alpha*(qhp/qcp+qhm/qcm - qhm/qcm*sqrt((Q-1).^2-4*qhp*mup*mum./(qhm*qcp^2*Beta.^2)))...
%        ./(2*(Beta+mup*mum./(qcm*qcp*Beta)))-d;   
mstar =@(alpha,Beta) pstar(alpha,Beta).*(qcm*Beta-qhm*alpha./(d+pstar(alpha,Beta)))/mum; 
Mstar = mstar(ALPHA,BETA);

%%% Unstable
pu_star =@(alpha,Beta) alpha*(qhp/qcp+qhm/qcm - qhm/qcm*sqrt((Q-1).^2-4*qhp*mup*mum./(qhm*qcp^2*Beta.^2)))...
       ./(2*(Beta+mup*mum./(qcm*qcp*Beta)))-d;
Pu_star = pu_star(ALPHA,BETA);
% Pstar2 = alpha*(qhp/qcp+qhm/qcm - qhm/qcm*sqrt((Q-1).^2-4*qhp*mup*mum./(qhm*qcp^2*Beta.^2)))...
%        ./(2*(Beta+mup*mum./(qcm*qcp*Beta)))-d;   
mu_star =@(alpha,Beta) pu_star(alpha,Beta).*(qcm*Beta-qhm*alpha./(d+pu_star(alpha,Beta)))/mum; 
Mu_star = mu_star(ALPHA,BETA);

Ju = J(Pu_star,Mu_star);
[Vu,Du] = eig(Ju);
epsilon = 1e-3;
%%% ODE solver stable manifold
Y0mp = [Pu_star,Mu_star]+epsilon*Vu(:,1)';
Y0mm = [Pu_star,Mu_star]-epsilon*Vu(:,1)';
Tfin = 1e3;
options = odeset('RelTol',1e-4,'AbsTol',1e-6);
%%% AMF induce competition on plant
[Tmp,Ymp] = ode45(@(t,y) -Func_1AMF_1Plant(y), 0:.1:Tfin, Y0mp, options);
[Tmm,Ymm] = ode45(@(t,y) -Func_1AMF_1Plant(y), 0:.1:Tfin, Y0mm, options);
UM = [[Ymp(end:-1:1,1);Pu_star;Ymm(:,1)],[Ymp(end:-1:1,2);Mu_star;Ymm(:,2)]];
%% Initial data
choice_init = 1;
if (choice_init == 1)
    %%%% Identical
    p0   = 2;
    m0   = .2;
elseif (choice_init == 2)
    %%% Close to (P*1,M*1)
    p0   = Pstar(1);  %.2;
    m0   = Mstar(1); %.2;
end
Y0 = [p0,m0];
%% ODE solver
Tfin = 1e3;
options = odeset('RelTol',1e-4,'AbsTol',1e-6);
%%% AMF induce competition on plant
[T,Y] = ode45(@(t,y) Func_1AMF_1Plant(y), 0:.1:Tfin, Y0, options);

%% Phase plane
D = qcm*BETA/(qhm*ALPHA);
ppmax = ALPHA*qhp./(qcp*BETA)-1;
mmmax = (Q-D)*(Q-1)/(mum*Q);
pp = 0:0.1:ppmax;
mm = 0:0.1:mmmax;
[PP,MM] = meshgrid(pp,mm);
dPP = Fp(PP,MM); 
dMM = Fm(PP,MM); 
nPM = sqrt(dPP.^2+dMM.^2);
dPP = .2*dPP./nPM;
dMM = .2*dMM./nPM;

%%%% Isoclines
p = 0:0.01:ppmax;
IPP =   mup*p./(ALPHA.*qhp./(d+p) - qcp*BETA);
IMM = (qcm*BETA*p - qhm*ALPHA.*p./(d+p))./mum;


figure(111)
clf
quiver(PP,MM,dPP,dMM)
hold on
plot(p,IPP,'linewidth',1.5,'color','b')
line([0,0],[0,mmmax],'color','b')
plot(p,IMM,'linewidth',1.5,'color','r')
line([0,ppmax],[0,0],'color','r')
plot(Y(:,1),Y(:,2),'Linewidth',1.5)
plot(UM(:,1),UM(:,2),'Linewidth',1.5,'color','k')

axis([0,ppmax,0,mmmax])



%% Figures
figure(1)
clf
hold on
plot(T,Y(:,1),'Linewidth',1.5)
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




