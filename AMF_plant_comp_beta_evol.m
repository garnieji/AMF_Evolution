%Write preliminary model for mycorrhiza
clear all
close all

options = odeset('RelTol',1e-4,'AbsTol',1e-6);

%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Time 
Tf = 2e2;
dt = 5e-3;
tt = 0:dt:Tf; Nt = length(tt); 

%%% Finctionnal parameters
% global rH qhp qcp mup mum alpha qcm qhm d ap12 ap21 disp_p disp_m

qhp = 3;
qcm = 2;
qcp = 1;
qhm = 1;
Q   = qcm*qhp/(qcp*qhm);

mup = 0.3;
mum = 0.3;

alpha = 5; %0.3; %1.5; 5

%%% Plant traits beta
betamax = 0.6; %.6 % .4
betamin = 0.4; % .4 % .4
dbeta   = 0.1;
Beta    = betamin:dbeta:betamax;
BETA    = Beta';
Nbeta   = length(Beta);

% Diffusion matrix
e = ones(Nbeta,1);
I_beta  = spdiags(e,0,Nbeta,Nbeta);
Ad_beta = spdiags([e -2*e e],-1:1,Nbeta,Nbeta);
Ad_beta(1,1) = -1;
Ad_beta(end,end) = -1;
Ad_beta = Ad_beta/(dbeta^2);

disp_p = 0.001; %0.0005;
disp_m = 0.001; %0.0005; % Mutation rates

A_pn = I_beta - dt*disp_p*Ad_beta;
A_mn = I_beta - dt*disp_m*Ad_beta;


rH    = 0.0;
d     = 1.2;
ap12  = .1*dbeta;
ap21  = .1*dbeta;

%% Steady states
% Ap_inv = [ap21,ap12];
% pstar =@(alpha,Beta) alpha*(qhp/qcp+qhm/qcm + qhm/qcm*sqrt((Q-1).^2-4*qhp*mup*mum./(qhm*qcp^2*Beta.^2)))...
%        ./(2*(Beta+mup*mum./(qcm*qcp*Beta)))-d;
% Pstar = pstar(alpha,Beta);
% mstar =@(alpha,Beta) pstar(alpha,Beta).*(qcm*Beta-qhm*alpha./(d+pstar(alpha,Beta)))/mum; 
% Mstar = mstar(alpha,Beta);
% 
% Dab = @(M) alpha*Ap_inv./(Ap_inv+M) - qcp*d*Beta/qhp;
% D = @(M)  qcm*Beta*d./(qhm*alpha*Ap_inv./(Ap_inv+M));

%% Functionnnal
fp = @(p,m,beta)  qhp*rH*p + qhp*alpha*m.*(ap12)./(ap12+(sum(p)-p)*dbeta).*p./(d+p) - ...
    qcp*beta.*m.*p - mup*p.^2;
fm = @(p,m,beta)  qcm*beta.*m.*p - qhm*alpha.*m.*(ap21)./(ap21+(sum(p)-p)*dbeta).*p./(d+p) - mum*m.^2;

% fp = @(p,m,beta)  qhp*rH*p + qhp*alpha*m.*(ap12)./(ap12+sum(p)*dbeta).*p./(d+p) - ...
%     qcp*beta.*m.*p - mup*p.^2;
% fm = @(p,m,beta)  qcm*beta.*m.*p - qhm*alpha.*m.*(ap21)./(ap21+sum(p)*dbeta).*p./(d+p) - mum*m.^2;


%% Initialization
P    = zeros(Nbeta,Nt);
M    = zeros(Nbeta,Nt);
PTOT = zeros (1,Nt);   
%% Initial data
P0  = 0.09*rand(Nbeta,1);
m0  = 1/Nbeta; %0.5;
M0  = m0.*(P0>0);
Ptot= sum(P0.*dbeta,1);

it = 1;
P(:,it) = P0; 
M(:,it) = M0;

Pnew = P0; Mnew = M0;
while (it<Nt)&&(Ptot>0)
    Pold = Pnew; Mold = Mnew;
    Pn = Pold + dt*fp(Pold,Mold,BETA);
    Mn = Mold + dt*fm(Pold,Mold,BETA);
    Pnew = A_pn\Pn;
    Mnew = A_mn\Mn;
    Ptot= sum(Pnew.*dbeta,1);

    it = it + 1;
    P(:,it) = Pnew; M(:,it) = Mnew;
    PTOT(it) = Ptot;
end


%% Figure
Color = get(gca,'colororder');
Marker = ['o','*','d','^'];
%%% No competition
% figure(1)
% for it = 1:1:Nt
%     clf
% 
%     hold on
%     yyaxis left
%     plot(Beta,P(:,it),'--')
%     
%     yyaxis right
%     plot(Beta,M(:,it),'--')
%     drawnow
%     hold off
%     pause(0.2)
% end




figure(1)
clf
hold on
for i = 1 :Nbeta
plot(tt,P(i,:),'--','color',Color(i,:))
% plot(tt,P(end,:),'--','color',Color(2,:))
plot(tt,M(i,:),'-','color',Color(i,:))
% plot(tt,M(end,:),'-','color',Color(2,:))
end


