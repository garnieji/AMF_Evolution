clear all
close all
% format LongE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Michorizer model with frequency dependence and evolution of alpha no space   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time parameter
TF = [70,100];
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
betamin = sqrt(4*Q*mup*mui./((Q-1).^2));
betamax = 13;
BETA = 0.72; %[0.43,0.72];
Nbeta = length(BETA);
% alpha = beta;
%% Trait alpha
alphamin = 0;
alphamax = 5;
dalpha = 0.01;                      % stepsize in alpha
ALPHA  = alphamin:dalpha:alphamax;  % vector of alpha values
Nalpha = length(ALPHA);
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
xmax = 40;
dx = 0.05;
xx = xmin:dx:xmax;
Nx = length(xx);

% Diffusion matrix
e = ones(Nx,1);
I_x  = spdiags(e,0,Nx,Nx);
Ad_x = spdiags([e -2*e e],-1:1,Nx,Nx);
Ad_x(1,1) = -1;
Ad_x(end,end) = -1;
Ad_x = Ad_x/(dx^2);

D_p = 0.1;
D_m = 0.1;  % diffusion rate

%% Spreading speed / Equilibrium
C_simu = zeros(2,Nbeta);
PSTAR_approx = zeros(1,Nbeta); MSTAR_approx = zeros(1,Nbeta);ABARSTAR_approx = zeros(1,Nbeta);
PSTAR = zeros(1,Nbeta); MSTAR = zeros(1,Nbeta);ABARSTAR = zeros(1,Nbeta);
P = zeros(Nx,Nbeta);
M = zeros(Nalpha,Nx,Nbeta);
parfor ib = 1:Nbeta
    beta  = max(betamin,BETA(ib));
%     beta = 8;
if (beta <2)
    Tf = TF(1);
else
    Tf = TF(2);
end
    %% Functional response and interactions MARIA
    fp = @(alpha,P,M)  P.*( q_hp*rp + (Q*sum(ALPHA'.*M*dalpha)./(d+P) -beta*sum(M*dalpha))-mup*P);
    fm = @(alpha,P,M)  M.*( (beta - alpha./(d+P)).*P - mui.*sum(M).*dalpha);
    
    %% Equilibrium
    Pstar  = (Q+1+sqrt( (Q-1)^2 -4*Q*mup*mui./(beta.^2)) )./(2*(beta+mup*mui./beta));
    z0  = -fzero(@(x) airy(1,x),0);
    PA = [1,-d/Pstar,0,-(z0)^3*dm];
    abar_approx = max(abs(roots(PA)));
    pstar    = abar_approx*Pstar-d;
    mstar = pstar.*(beta-1./Pstar)/mui;
    eta = (pstar./(pstar+d)/dm)^(1/3);
    m_d_star = airy( eta*ALPHA'-z0 );
    Mstar = m_d_star./sum(m_d_star*dalpha)*mstar;  
    PSTAR_approx(ib) = pstar; MSTAR_approx(ib) = mstar; ABARSTAR_approx(ib) = abar_approx;
    %% Initial data
    P0 = pstar*(xx<=0);
    M0 = Mstar.*(xx<=0);
    
    t = 0; it = 1;
    dt = 0.01;
    tt = 0:dt:Tf;
    Nt = length(tt);
    Pnew = P0;  Mnew = M0;
    
    Xt = zeros(2,Nt);
    while (it<Nt)
        Pold = Pnew;Mold = Mnew;
        
        Pp = Pold + dt*fp(ALPHA',Pold,Mold);
        Bp = I_x-dt*D_p*Ad_x;
        Ppnew = Bp\(Pp)';
        Pnew = Ppnew';
        
        Mm = Mold + dt*fm(ALPHA',Pold,Mold);
        Am = I_alpha-dt*dm*Ad_alpha;
        Mnew = Am\(dt*D_m*Mold*Ad_x+Mm);
        
       
        it = it + 1; t = t+dt;
        xpt = max(sum(Pnew>(pstar/2)),1);
        xmt = max(sum(sum(Mnew*dalpha,1)>(mstar/2)),1);
        Xt(1,it) = xx(xpt);
        Xt(2,it) = xx(xmt);
    end
        Abar = ALPHA*Mnew./sum(Mnew,1);
        PSTAR(ib) = mean(Pnew(1:100));
        MSTAR(ib) = mean(sum(Mnew(:,1:100),1)*dalpha); 
        ABARSTAR(ib) = mean(Abar(1:100));
%     pxt = polyfit(tt(end-500:end),Xt(end-500:end),1);
%     c_simu = pxt(1);
    c_simu = mean(Xt(:,end-500:end)/(tt(end-500:end)),2);
    C_simu(:,ib) = c_simu;
    P(:,ib) = Pnew;
    M(:,:,ib) = Mnew;
end
M_x = permute(sum(M*dalpha,1),[2,3,1]);
Abar = permute(sum(ALPHA'.*M,1)./sum(M,1),[2,3,1]);
%% Figure
Color = get(gca,'colororder');

bluegradient = ["#00008b"; "#1d289b"; "#314fac"; "#5877bd"; "#759fcd"; "#92c6de"; "#afeeee"];

bluegradient2 = [ [97 148 188]/255; [116 162 199]/255; [134 177 210]/255; [152,191,222]/255; ...
    [171 205 233]/255; [190 220 244]/255; [208 234 255]/255];
cmp = colormap('winter');

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



bluegradient3 = [flipud([cmp(1:5:end-50,:);cmp(end-40:20:end,:) ]);...
    [97 148 188]/255; [116 162 199]/255; [134 177 210]/255; [152,191,222]/255; ...
    [171 205 233]/255; [190 220 244]/255; [208 234 255]/255];


red2 = [188 70 118]/255;
redgradient = [ red2; [198 98 139]/255; [208 126 160]/255; [218 136 170]/255];

%% Figure Front 

aa = [0,0.25,0.5,0.75,1,2,3];
x = [0,10,15,30];
% x = [0,6,12,30];
Ixx = sum(xx<=x',2)';
Nxx = length(x);
Iaa = sum(ALPHA<=aa',2)';
Naa = length(aa);
ib = 1;
Mc = cumsum(M*dalpha,1);
p = zeros(1,Naa);
figure(1)
clf
hold on
% yyaxis left
for ia=1:Naa
p(ia) = plot(xx,M(Iaa(ia),:,ib),'LineWidth',2,'Color',outpict(ia,:),...
    'LineStyle','-','Marker','none','DisplayName',['$\alpha=$',num2str(ALPHA(Iaa(ia)))]);
end
% for ia=1:Naa
% p(ia) = plot(xx,Mc(Iaa(ia),:,ib),'LineWidth',2,'Color',outpict(ia,:),...
%     'LineStyle','-','Marker','none','DisplayName',['$\alpha=$',num2str(ALPHA(Iaa(ia)))]);
% end
plot(xx,P(:,ib),'Color','k','LineWidth',2,'LineStyle','--')
plot(xx,M_x(:,ib),'Color','k','LineWidth',2,'LineStyle','-')

ylabel(['Densities of ','AMF and plant'],'Interpreter','latex','FontSize',20)


yyaxis right
for ib=1:Nbeta
plot(xx,Abar(:,ib),'Color',red2,'LineWidth',2)
% Ax = z0*(dm*(P(:,ib)+d)./P(:,ib)).^(1/3);
% plot(xx,Ax,'--')
% line([xx(1),xx(end)],[ABARSTAR_approx(ib),ABARSTAR_approx(ib)],'linestyle','--','color',red2)
% line([xx(1),xx(end)],[d/Pstar,d/Pstar],'linestyle','--','color',red2)
end
for ix = 1:Nxx
line([x(ix),x(ix)],[0,Abar(Ixx(ix),ib)],'linestyle',':','color',[0.85,0.85,0.85],'linewidth',1.5)
plot(x(ix),Abar(Ixx(ix),ib),'+k', 'MarkerSize', 7,'LineWidth',2); 
end
ylim([0,1])
xlim([xmin,40])
% hl = legend(p,'location','southeast');
% set(hl, 'Interpreter','latex', 'FontSize',16)
ax = gca;
c = ax.YColor;
ax.YColor = red2;
ylabel({'Mean mutualistic investment', 'of AMF community, $\overline{\alpha}$'},'Interpreter','latex','FontSize',20)
xlabel('space, $x$','Interpreter','latex','FontSize',20)

%% Alpha mean
M_d = M./sum(M*dalpha);

redgradient = [ red2; [198 98 139]/255; [208 126 160]/255; [218 136 170]/255];
CT = cbrewer('seq','PuRd',5);
redgradient = CT(2:end,:    )
style = {'-','--','-.',':'};
figure(2)
clf
hold on
p = zeros(1,Nxx);
for ix = 1:Nxx
    p(ix) = plot(ALPHA,M_d(:,Ixx(ix),ib),'Color',redgradient(ix,:),'LineWidth',2,...
                 'LineStyle',style{ix},'DisplayName',['$x=$',num2str(xx(Ixx(ix)))])
end
hl = legend(p);
set(hl, 'Interpreter','latex', 'FontSize',16)
ylabel({'Distribution of mutualistic','investment in the AMF community'},'Interpreter','latex','FontSize',20)
xlabel('AMF trait, $\alpha$','Interpreter','latex','FontSize',20)

%% Growth rate 
beta = BETA;
GP = ((Q*Abar./(P+d)-beta).*M_x-mup*P);
GM_i = ( (beta -ALPHA./(d+P)).*P - mui.*M_x); 
GM =  ( (beta -Abar./(d+P)).*P - mui.*M_x);


figure(20)
clf
line([xmin,xmax],[0,0],'color','k')
hold on
plot(xx,P,'--k','LineWidth',2)
plot(xx,M_x,'-k','LineWidth',2)
ylabel(['Densities of ',' AMF and plant'],'Interpreter','latex','FontSize',20)
% axis([-5,20,-0.005,1.4])

yyaxis right
plot(xx,GM,'-','LineWidth',2,'color',red2)
hold on
plot(xx,GP,'--','LineWidth',2,'color',red2)
ylabel({'Growth rate of','plant and AMF community'},'Interpreter','latex','fontsize',20)
axis([0,20,0,0.12])
% axis([-5,20,-0.005,0.06])
% ylabel({'Growth rate of','Plant and total AMF'},'Interpreter','latex','fontsize',20)
xlabel('space, $x$','Interpreter','latex','fontsize',20)

ax = gca;
ax.YAxis(2).Color = red2;
% ax.YAxis(1).Limits(1) = ax.YAxis(2).Limits(2); 

%%
M_d = M./sum(M*dalpha);
figure(3)
clf
hold on
p = zeros(1,Nbeta);
for ib = 1:Nbeta
    p(ib) = scatter(ALPHA,M_d(:,1,ib),'Filled','MarkerFaceColor',redgradient(ib,:),'DisplayName',['\beta=',num2str(BETA(ib))])
    scatter(ALPHA,M_d(:,201,ib),'d','Filled','MarkerFaceColor',redgradient(ib,:))
    scatter(ALPHA,M_d(:,401,ib),'v','Filled','MarkerFaceColor',redgradient(ib,:))
end
legend(p,'FontSize',16)

%% Parasitic vs Mutualistic
Gr = cbrewer('seq','Greens',8);
beta = BETA;
Pstar  = (Q+1+sqrt( (Q-1)^2 -4*Q*mup*mui./(beta.^2)) )./(2*(beta+mup*mui./beta));
alpha_c = max(beta*d/Q,d/Pstar);
Im_parasitic   = (ALPHA<=alpha_c);
Im_mutualistic = (ALPHA>alpha_c);
M_p = M(Im_parasitic,:);
Mx_p = sum(M_p*dalpha);
Mx_m =sum(M(Im_mutualistic,:)*dalpha);
Prop_parasitic = sum(M_p)./sum(M)*100;

figure(4)
clf
hold on
xxx = [xx, fliplr(xx)];
inBetween = [zeros(1,length(xx)), fliplr(Mx_p)];
f = fill(xxx, inBetween,outpict(1,:));%red2);
set(f,'EdgeColor','none','FaceAlpha', 0.5)
plot(xx,Mx_p,'LineWidth',2,'color',outpict(1,:))%red2)
ylim([0,1.8])

text(-4.5,0.7,'Mutualistic','Interpreter','latex','Color','b','FontSize',17,'FontWeight','bold')


inBetween = [Mx_p, fliplr(M_x')];
f = fill(xxx, inBetween,outpict(end,:));
set(f,'EdgeColor','none','FaceAlpha', 0.5)
plot(xx,M_x,'k','LineWidth',3)

text(-3.8,0.4,'Parasitic','Interpreter','latex','Color',Gr(8,:),'FontSize',20)

ylabel({'Cumulative densities of parasitic','and mutualistic AMF'},'Interpreter','latex','FontSize',20)

yyaxis right
plot(xx,Prop_parasitic,'color',red2,'LineWidth',2)
ax = gca;
c = ax.YColor;
ax.YColor = red2;
ylim([0,40])
% plot(xx,Mx_m+Mx_p,'--')
ylabel({'Percentage of parasitic AMF','along the front'},'Interpreter','latex','FontSize',20)
xlabel('space, $x$','Interpreter','latex','FontSize',20)



