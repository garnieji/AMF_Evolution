clear all
close all
% format LongE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Michorizer model with frequency dependence                           %
%        evolution of alpha spreading in space                            %
%          LINEAGES                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Parameter of the model
% global q_hp q_cm q_hm q_cp beta mup mui d rp
q_hp  = 3; % q>1
q_cm  = 2;
q_hm  = 1;
q_cp  = 1;
Q = q_hp*q_cm/(q_cp*q_hm);

mup = 0.3; % 1/100
mui = 0.3; % 1/20

d = 1.2;                            % plant saturation parameter
%% Plant
rp = 0; %.02;
% betamin = sqrt(4*Q*mup*mui./((Q-1).^2));
% betamax = 13;
% BETA = [linspace(betamin,2,25),linspace(2.1,betamax,15)];
beta = 0.7;
% Nbeta = length(BETA);
% alpha = beta;
%% Trait alpha
alphamin = 0;
alphamax = 5;
dalpha = 0.01;                      % stepsize in alpha
ALPHA  = alphamin:dalpha:alphamax;  % vector of alpha values
Nalpha = length(ALPHA);
%% Space x
xmin = -5;
xmax = 20;
dx = 0.05;
xx = xmin:dx:xmax;
Nx = length(xx);
%% Lineages
yy = linspace(xmin,5,10);
Ny = length(yy);
aa = linspace(alphamin,alphamax,10);
Na = length(aa);

%% Load Data
% load('Evol_lineages.mat')
% load('Evol_lineages_1.mat')
load('Evol_lineages2_iy_1.mat')
Nx = length(xx);
Ny = length(yy);
Na = length(aa);

v = zeros(Nalpha,Nx,Na,Ny);
vinit = zeros(Nalpha,Nx,Na,Ny);
vp = zeros(1,Nx,Ny);
vpinit = zeros(1,Nx,Ny);
 
for iy = 1:Ny
    load(['Evol_lineages2_iy_',num2str(iy),'.mat'])
    v(:,:,:,iy) = V;
    vinit(:,:,:,iy) = Vinit;
    vp(:,:,iy) = Vp;
    vpinit(:,:,iy) = Vpinit;
end

%% Mnew 
V_x = permute(sum(v*dalpha,1),[4,3,2,1]);
V_x = reshape(V_x,Na*Ny,Nx);
Vx = cumsum(V_x,1);
Vp_x = permute(vp,[3,2,1]);
M = sum(v,[3,4]);
Abar = sum(ALPHA'.*M)./sum(M);
M_x = sum(M*dalpha);


% figure(1)
% clf
% plot(xx,M_x,'--','LineWidth',5)
% hold on
% % plot(xx,sum(Minit*dalpha))
% % if (iy<Ny)
% % plot(xx,sum(Minit.*((xx<yy(iy+1)).*(xx>=yy(iy)))*dalpha))
% % else
% % plot(xx,sum(Minit.*((xx>yy(iy)))*dalpha))
% % end 
% for ii = 1:Ny*Na
%         plot(xx,Vx(ii,:),'-')
% end
%% Figure Proportion
Prop = V_x./M_x;
pr = reshape(Prop,Ny,Na,Nx);
for iy = 1:Ny
    for ia = 1:Na
         if (iy<Ny)&&(ia<Na)
            di = (yy(iy+1)-yy(iy)).*(aa(ia+1)-aa(ia));
        elseif(iy==Ny)&&(ia<Na)
            di = (xx(end)-yy(iy)).*(aa(ia+1)-aa(ia));
        elseif(iy<Ny)&&(ia==Na)
            di = (yy(iy+1)-yy(iy)).*(ALPHA(end)-aa(ia));
        elseif (iy==Ny)&&(ia==Na)
            di = (xx(end)-yy(iy)).*(ALPHA(end)-aa(ia));
        end
        pr(iy,ia,:) = pr(iy,ia,:)./(di+(di<=0)).*(di>0);
    end
end
Prop = reshape(pr,Na*Ny,Nx);


prop= cumsum(Prop);
phi = reshape(Prop(:,300),Ny,Na);
daa = aa(2)-aa(1);
dyy = yy(2)-yy(1);
phi = phi./trapz(aa,trapz(yy,phi,2));
% phi = phi./(sum(phi,'all')*daa*dyy);
% figure(2)
% clf
% plot(xx,prop)
% hold on
% plot(xx,P,'-k','linewidth',2)
% plot(xx,M_x,'--k','linewidth',2)

%% Growth rate
beta = 0.7;
Abarinit = sum(ALPHA'.*Minit)./sum(Minit);
Minit_x = sum(Minit*dalpha);
GP = ((Q*Abarinit./(Pinit+d)-beta).*Minit_x-mup*Pinit);
GM_i = ( (beta -ALPHA'./(d+Pinit)).*Pinit - mui.*Minit_x); 
GM =  ( (beta -Abarinit./(d+Pinit)).*Pinit - mui.*Minit_x);

figure(20)
clf
yyaxis right
plot(xx,GM,'-','LineWidth',2)
hold on
plot(xx,GP,'--','LineWidth',2)
ylabel({'Growth rate of','Plant and total AMF'},'Interpreter','latex','fontsize',20)

yyaxis left
plot(xx,Pinit,'--','LineWidth',2)
hold on
plot(xx,Minit_x,'-','LineWidth',2)
xlim([yy(1),10])
ylabel('Density of Plant and AMF','Interpreter','latex','fontsize',20)
xlabel('space, $x$','Interpreter','latex','fontsize',20)

figure(2000)
clf
contourf(xx,ALPHA,GM_i)
hold on
yyaxis right
plot(xx,Pinit,'--','LineWidth',2)
hold on
plot(xx,Minit_x,'-','LineWidth',2)
yyaxis left
plot(xx,GM,'-','LineWidth',2)
hold on
plot(xx,GP,'--','LineWidth',2)

figure(200)
clf
scatter(xx(1:5:end),Pinit(1:5:end),15,GP(1:5:end),'filled','^')
hold on
scatter(xx,Minit_x,[15],GM,'filled')
xlim([yy(1),10])
ylabel('Density of Plant and AMF','Interpreter','latex','fontsize',20)
xlabel('space ($x$)','Interpreter','latex','fontsize',20)
colormap winter



ym = linspace(0,max(Minit_x)+0.1,10);
yp = linspace(0,max(Pinit)+0.1,10);
GMM = ones(length(ym),1)*GM;
GPP = ones(length(yp),1)*GP;
figure(201)
clf
hold on
contourf(xx,ym,GMM,200,'edgecolor','none')
plot(xx,Minit_x,'k--','LineWidth',2)
colormap winter
xlim([yy(1),10])

figure(202)
clf
hold on
contourf(xx,yp,GPP,200,'edgecolor','none')
plot(xx,Pinit,'k','LineWidth',4)
colormap winter
xlim([yy(1),10])
ylabel('Population density','Interpreter','latex','fontsize',20)
xlabel('position ($x$)','Interpreter','latex','fontsize',20)
%% 
mSTAR = M(:,10);
MSTAR = sum(mSTAR)*dalpha;

% Minit = zeros(Nalpha,Nx); 
% NNx = sum(sum(M*dalpha)<MSTAR);
% Nxx = Nx - NNx+1;
% Minit(:,1:NNx) = M(:,Nxx:end);
Minit = Minit./sum(Minit*dalpha*dx,'all');
ivx = sum(M_x>M_x(1)/2);
VVx = reshape(V_x(:,ivx),Ny,Na);

Vvpx = Vp_x(:,ivx);

Prop_v = VVx./M_x(ivx);
Prop_V = zeros(Ny,Na);
for iy = 1:Ny
    for ia = 1:Na
         if (iy<Ny)&&(ia<Na)
            di = (yy(iy+1)-yy(iy)).*(aa(ia+1)-aa(ia));
        elseif(iy==Ny)&&(ia<Na)
            di = (xx(end)-yy(iy)).*(aa(ia+1)-aa(ia));
        elseif(iy<Ny)&&(ia==Na)
            di = (yy(iy+1)-yy(iy)).*(ALPHA(end)-aa(ia));
        elseif (iy==Ny)&&(ia==Na)
            di = (xx(end)-yy(iy)).*(ALPHA(end)-aa(ia));
        end
        Prop_V(iy,ia) = Prop_v(iy,ia)./(di+(di<=0)).*(di>0);        
    end   
end

IVVx = trapz(yy,trapz(aa,VVx,2));
VVx = VVx./IVVx;
IVVpx = trapz(yy,Vvpx);
VVpx = Vvpx./IVVpx;
Abar_int = sum(ALPHA'.*Minit)./sum(Minit);

figure(21)
clf
contourf(yy,aa,VVx')
hold on
contour(xx,ALPHA,Minit,'LineWidth',2)
xlabel('space x')
ylabel('trait \alpha')
% contour(xx,ALPHA,Minit.^2./sum(Minit.^2*dalpha*dx,'all'),'LineWidth',2,'LineStyle','--')

xlim([yy(1),yy(end)])


% [XX,ALpha]= meshgrid(xx,ALPHA);
% [XXa,AAx]= meshgrid(xx,aa);
% [YY,AA] = meshgrid(yy,aa);
% VVxalpha = interp2(YY,AA,VVx',XXa,AAx);
% Minit_ya = interp2(XX,ALpha,Minit,XXa,AAx);
% 
% figure(210)
% clf
% hold on
% for ia = 1:5
% scatter(xx,Minit_ya(ia,:),[],VVxalpha(ia,:),'filled')
% end
% colormap winter
% % contour(xx,ALPHA,Minit.^2./sum(Minit.^2*dalpha*dx,'all'),'LineWidth',2,'LineStyle','--')
% 
% xlim([yy(1),yy(end)])

VVy = trapz(aa,VVx,2);
VVay = trapz(aa,aa.*VVx,2)./VVy;
ix = sum(xx<=yy(end));
VVy = interp1(yy,VVy,xx(1:ix),'makima');
VVay = interp1(yy,VVay,xx(1:ix),'makima');
PVy = VVy./sum(VVy);
ix = sum(xx<=yy(end));
VVpy = interp1(yy,VVpx,xx(1:ix),'makima');

cm = [ones(ix,1) xx(1:ix)']\log(VVy./(Minit_x(1:ix).^2))';
cp = [ones(ix,1) xx(1:ix)']\log(VVpy./(Pinit(1:ix).^2))';
c= mean([cm(2),cp(2)]);
VVy_approx = exp(c.*xx(1:ix)).*(Minit_x(1:ix).^2)./sum(exp(c.*xx(1:ix)).*(Minit_x(1:ix).^2)*dx);


VVpy_approx = exp(c.*xx(1:ix)).*(Pinit(1:ix).^2)./sum(exp(c.*xx(1:ix)).*(Pinit(1:ix).^2)*dx);
% VVpy_approx = exp(C_simu(7).*xx(1:ix)).*(Pinit(1:ix).^2)./sum(exp(C_simu(7).*xx(1:ix)).*(Pinit(1:ix).^2)*dx);

% Minit_x = sum(Minit*dalpha);
%%
% basic 2-point interpolation
colorA = [0.1,0.8,0.6];
colorB = [97 148 188]/255;
sz = [3 3]; % [y x]
x0 = 1/sz(1);
xq = linspace(x0,1,sz(1));
outpict = interp1([x0 1],[colorA; colorB],xq,'linear','extrap');
colorA = [116 162 199]/255;
colorB = [0.1 0.25 .8];
sz = [4 4]; % [y x]
x0 = 1/sz(1);
xq = linspace(x0,1,sz(1));
outpict = [outpict;interp1([x0 1],[colorA; colorB],xq,'linear','extrap')];
%%
figure(22)
clf
yyaxis right
plot(xx(1:ix),VVy_approx,'-','LineWidth',3,'color',[0.8,0.8,0.8])
hold on
% plot(xx(1:ix),VVy,'-','LineWidth',2)
scatter(xx(1:3:ix),VVy(1:3:end),25,VVay(1:3:end),'d','filled')

% scatter(xx(1:ix),VVy,100,Abarinit(1:ix),'filled')
colormap(outpict)


plot(xx(1:ix),VVpy_approx,'--','LineWidth',3,'color',[0.8,0.8,0.8])
plot(xx(1:ix),VVpy,'--','LineWidth',2)
% plot(xx(1:ix),Abarinit(1:ix),'LineWidth',2)
% plot(xx,GM,'-','LineWidth',2)
% plot(xx,GP,'--','LineWidth',2)
% plot(xx(1:ix),PVy,'LineWidth',2)
% plot(yy,VVy,'','LineWidth',2)

ylabel({'Fixation probability of AMF','and plant common ancestors'},'Interpreter','latex','fontsize',20)
hold on
yyaxis left
plot(xx,Pinit,'k--','LineWidth',2)
% plot(xx,Minit_x,'k','LineWidth',2)
scatter(xx,Minit_x,[20],Abarinit,'filled')

% scatter(xx,Pinit,[15],GP,'filled','Marker','d')
% colormap('gray')
xlim([yy(1),8])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
caxis([0,0.8])
colorbar('Location','northoutside')
ylabel('Densities of AMF and plant','Interpreter','latex','fontsize',20)
xlabel('space, $x$','Interpreter','latex','fontsize',20)

%%

% yyaxis right
% plot(xx,GM,'-','LineWidth',2)
% hold on
% plot(xx,GP,'--','LineWidth',2)
% ylabel({'Growth rate of','Plant and total AMF'},'Interpreter','latex','fontsize',20)


%%

VVa =  VVx(16,:)./trapz(aa,VVx(16,:));
VVa = interp1(aa,VVa,ALPHA,'makima');

Minit_a = Minit(:,1)./sum(Minit(:,1)*dalpha);
figure(23)
clf
% plot(aa,VVa)
scatter(ALPHA,VVa,'filled')

hold on
plot(ALPHA,Minit_a,'LineWidth',2)
plot(ALPHA,Minit_a.^2./sum(Minit_a.^2*dalpha),'--','LineWidth',2)

%% Figure x,\alpha
figure(3)
clf
contour(xx,ALPHA,V(:,:,3,1)./M,'LineStyle','--')
hold on
contour(xx,ALPHA,V(:,:,6,1)./M)








