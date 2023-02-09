clear all
close all
% format LongE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Michorizer model with frequency dependence and evolution of alpha no space   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time parameter
Tf = 7e1;
%% Parameter of the model
global q_hp q_cm q_hm q_cp beta mup mui d rp Aa a
q_hp  = 3; % q>1
q_cm  = 2;
q_hm  = 1;
q_cp  = 1;

beta  = 0.4;
% alpha = beta;
%% Trait alpha
alphamin = 0;
alphamax = 1;
dalpha = 0.05;
ALPHA  = alphamin:dalpha:alphamax;
Nalpha = length(ALPHA);

mup = 0.3; % 1/100
mui = 0.3; % 1/20

d = 1.2;

rp = 0;%.02;
%% Choice case 3
    a = 0.2;
    a_ww = a;
    a_cw = a;
    a_wc = a;
    a_cc = a;

N_AMF = Nalpha;

% Diffusion matrix
e = ones(Nalpha,1);
I_alpha  = spdiags(e,0,Nalpha,Nalpha);
Ad_alpha = spdiags([e -2*e e],-1:1,Nalpha,Nalpha);
Ad_alpha(1,1) = -1;
Ad_alpha(end,end) = -1;
Ad_alpha = Ad_alpha/(dalpha^2);

dm = 0.01;  % mutation rate

%% Space x
xmin = -5;
xmax = 30;
dx = 0.05;
xx = xmin:dx:xmax;
Nx = length(xx);

% Diffusion matrix
e = ones(Nx,1);
I_x  = spdiags(e,0,Nx,Nx);
Ad_x = spdiags([e -2*e e],-1:1,Nx,Nx);
Ad_x(1,1) = -1;
% Ad_x(end,end) = -1;
Ad_x = Ad_x/(dx^2);

D_p = .1;
D_m = .1;  % diffusion rate

choice_comp = 3;
if (choice_comp == 1)
    %% Functionnal response without competition
    fp = @(alpha,P,M)  P.*( q_hp*rp + (q_hp*sum(alpha.*M*dalpha)./(d+P) -q_cp*beta*sum(M*dalpha) )-mup*P);
    fm = @(alpha,P,M)  M.*( (q_cm*beta - q_hm*alpha./(d+P)).*P - mui.*M);
elseif(choice_comp == 2)
    %% Functionnal response with competition on uptake (weak)
    fp = @(alpha,P,M)  P.*( q_hp*rp + (q_hp*sum(alpha.*M*dalpha)./(d+P) -q_cp*beta*a*sum(M*dalpha)./(a+sum(M*dalpha)) )-mup*P);
    fm = @(alpha,P,M)  M.*( (q_cm*beta*a./(a+sum(M*dalpha)) - q_hm*alpha./(d+P)).*P - mui.*M);
elseif (choice_comp == 3)
    %% Functionnal response with competition on maintenance cost (strong)
    fp = @(alpha,P,M)  P.*( q_hp*rp + (q_hp*sum(alpha.*M*dalpha)./(d+P) -q_cp*beta*sum(M*dalpha))-mup*P);
    fm = @(alpha,P,M)  M.*( (q_cm*beta - q_hm*alpha./(d+P)).*P - mui.*sum(M).*dalpha);
elseif (choice_comp == 4)
    %% Functionnal response with competition on maintenance cost and uptake (weak+strong)
    fp = @(alpha,P,M)  P.*( q_hp*rp + (q_hp*sum(alpha.*M*dalpha)./(d+P) -q_cp*beta*a*sum(M*dalpha)./(a+sum(M*dalpha)) )-mup*P);
    fm = @(alpha,P,M)  M.*( (q_cm*beta*a./(a+sum(M*dalpha)) - q_hm*alpha./(d+P)).*P - mui.*sum(M).*dalpha);
end
F = @(X) [fp(ALPHA',X(1),X(2:end));...
          dm*Ad_alpha*X(2:end)+fm(ALPHA',X(1),X(2:end))];
%% Initial data
p0 = 9;
m0 = 4*ones(1,N_AMF);
Tfe = 1e2;
[T,X] = ode45(@(t,X) F(X),[0,Tfe],[p0,m0]);
Pinfty = X(end,1);
Minfty = X(end,2:end);
P0 = Pinfty*(xx<=0);
M0 = Minfty'.*(xx<=0);% ((ALPHA'<=0.1).*(ALPHA'>=0))

t = 0; it = 1; 
dt = 0.01;
tt = t:dt:Tf;
Nt = length(tt);
%% Variables
PP = zeros(Nt,Nx); MM = zeros(Nalpha,Nx,Nt);
MM_x = zeros(Nt,Nx); MM_d = zeros(Nt,Nalpha); MM_b = zeros(1,Nt);
Pnew = P0;  PP = P0;  
Mnew = M0;
PP(1,:) = P0;
MM(:,:,1) = M0;
MM_x(1,:) = sum(M0)*dalpha; MM_d(1,:) = (sum(M0,2)*dx)';MM_b(1) = sum(sum(M0)*dalpha)*dx;
for it = 2:Nt
    Pold = Pnew;Mold = Mnew;
    
    Pp = Pold + dt*fp(ALPHA',Pold,Mold);
    Bp = I_x-dt*D_p*Ad_x;
    Pnew = Pp/Bp;
    
%     Mm = Mold + dt*fm(ALPHA',Pold,Mold);
%     Am = I_alpha-dt*dm*Ad_alpha;
%     Bm = I_x-dt*D_m*Ad_x;
%     Mnew = Mold/Bm;
%     Mnew=Am\Mnew;
    
    Mm = Mold + dt*fm(ALPHA',Pold,Mold);
    Am = I_alpha-dt*dm*Ad_alpha;
    Bm = I_x-dt*D_m*Ad_x;
%     Mnew = Am\(dt*D_m*Mold*Ad_x+Mm);
    Mnew = (dt*dm*Ad_alpha*Mold + Mm)/Bm;
 
    PP(it,:) = Pnew;
    MM(:,:,it) = Mnew;
    MM_xt = sum(Mnew)*dalpha;
    MM_dt = sum(Mnew,2)*dx;
    MM_d(it,:) = MM_dt';
    MM_x(it,:) = MM_xt;
    MM_b(it) = sum(MM_xt)*dx;
end

%% Plot biommass 
PP_b = sum(PP,2)*dx;

figure(1)
clf
hold on
plot(tt,PP_b,'--')
plot(tt,MM_b,'-o')
xlim([0,Tf])
drawnow
hold off

%% Plot of M over space
% figure(2)
% clf
% for It = 1:50:length(tt)
%     figure(2)
%     clf
%     hold on
%     plot(xx,MM_x(It,:),'-')
%     plot(xx,PP(It,:),'--')
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

%% Plot of M distribution over time over space trait
MM_d = MM_d./MM_b';
mean_alpha = permute(sum(ALPHA'.*MM,1)./sum(MM,1),[2,3,1]);
% figure(2)
% clf
% plot(aalpha,MM_d(end,:))
It = Nt;
% for It = 1:50:length(tt)
% %     mean_alpha = sum(ALPHA'.*MM(:,:,It)*dalpha)./sum(MM(:,:,It)*dalpha);
%     figure(3)
%     clf
%     hold on
% %     plot(xx,mean_alpha)
%     plot(ALPHA,MM(:,10,It)./sum(MM(:,10,It)*dalpha),'-')
%     plot(ALPHA,MM(:,301,It)./sum(MM(:,301,It)*dalpha),'-')
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

figure(2)
clf
hold on
plot(ALPHA,MM(:,10,end)./sum(MM(:,10,It)*dalpha),'-')
plot(ALPHA,MM(:,301,end)./sum(MM(:,301,It)*dalpha),'--')
ylabel('Trait distribution of AMF $m(t,x,\alpha)$','Interpreter', 'latex','FontSize',16)
xlabel('trait $\alpha$','Interpreter','latex','FontSize',16)

%%
for It = 1:50:length(tt)
figure(3)
clf
hold on
yyaxis left
plot(xx,mean_alpha(:,It))
ylabel('Mean trait of AMF $\displaystyle mean(\alpha)(x) =\int_{\overline\alpha}^{\underline\alpha} {\alpha\,m(t,x,\alpha)\,d\alpha}$','Interpreter', 'latex','FontSize',16)
yyaxis right
plot(xx,MM_x(It,:),'-')
% plot(xx,MM(2,:,It),'-o')
% plot(xx,MM(20,:,It),'-d')
plot(xx,PP(It,:),'--')
ylabel('Density of AMF and plant ','Interpreter', 'latex','FontSize',16)
xlabel('space $x$','Interpreter','latex','FontSize',16)
end


