clear all
close all
% format LongE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Michorizer model with frequency dependence and evolution of alpha no space   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time parameter
Tf = 70;
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
xmax = 20;
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

%% Initialization Equilibrium
%% Functional response and interactions MARIA
fp = @(alpha,P,M)  P.*( q_hp*rp + (Q*sum(ALPHA'.*M*dalpha)./(d+P) -beta*sum(M*dalpha))-mup*P);
fm = @(alpha,P,M)  M.*( (beta - alpha./(d+P)).*P - mui.*sum(M).*dalpha);
fv = @(alpha,P,M)  ( (beta - alpha./(d+P)).*P - mui.*sum(M).*dalpha);
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
%% Initial data
P0 = pstar*(xx<=0);
M0 = Mstar.*(xx<=0);

t = 0; it = 1;
dt = 0.01;
tt = 0:dt:Tf;
Nt = length(tt);
Pnew = P0;  Mnew = M0;

while (it<Nt)
    Pold = Pnew;Mold = Mnew;
    
    Pp = Pold + dt*fp(ALPHA',Pold,Mold);
    Bp = I_x-dt*D_p*Ad_x;
    Ppnew = Bp\(Pp)';
    Pnew = Ppnew';
    
    Mm = Mold + dt*fm(ALPHA',Pold,Mold);
    Am = I_alpha-dt*dm*Ad_alpha;
    Mnew = Am\(dt*D_m*Mold*Ad_x+Mm);
    
    it = it + 1;
end
PSTAR = Pnew(10);
mSTAR = Mnew(:,10);
MSTAR = sum(mSTAR)*dalpha;
P = Pnew;
M = Mnew;

Pinit = zeros(1,Nx); Minit = zeros(Nalpha,Nx); 
NNx = sum(sum(M*dalpha)<MSTAR);
Nxx = Nx - NNx+1;
Pinit(1:NNx) = P(Nxx:end);
Minit(:,1:NNx) = M(:,Nxx:end);

%% Lineages
yy = linspace(xmin,5,10);
Ny = length(yy);
aa = linspace(alphamin,alphamax,10);
Na = length(aa);
Tf = 80;
tt = 0:dt:Tf;
Nt = length(tt);
V = zeros(Nalpha,Nx,Ny,Na);
for iy = 1:Ny
    parfor ia = 1:Na
        if (iy<Ny)&&(ia<Na)
            vi = Minit.*((xx<yy(iy+1)).*(xx>=yy(iy))).*(ALPHA'<aa(ia+1)).*(ALPHA'>=aa(ia));
        elseif(iy==Ny)&&(ia<Na)
            vi = Minit.*(xx>yy(iy)).*(ALPHA'<aa(ia+1)).*(ALPHA'>=aa(ia));
        elseif(iy<Ny)&&(ia==Na)
            vi = Minit.*((xx<yy(iy+1)).*(xx>=yy(iy))).*(ALPHA'>aa(ia));
        elseif (iy==Ny)&&(ia==Na)
            vi = Minit.*(xx>yy(iy)).*(ALPHA'>aa(ia));
        end
        it = 1;
        Pnew = Pinit;  Mnew = Minit; vnew = vi;
        while (it<Nt)
            Pold = Pnew;Mold = Mnew;vold = vnew;
            
            Pp = Pold + dt*fp(ALPHA',Pold,Mold);
            Bp = I_x-dt*D_p*Ad_x;
            Ppnew = Bp\(Pp)';
            Pnew = Ppnew';
            
            Mm = Mold + dt*fm(ALPHA',Pold,Mold);
            Am = I_alpha-dt*dm*Ad_alpha;
            Mnew = Am\(dt*D_m*Mold*Ad_x+Mm);
            
            vm = vold + dt*vold.*fv(ALPHA',Pold,Mold);
            vnew = Am\(dt*D_m*vold*Ad_x+vm);
            
            it = it + 1;
%             if (mod(it,10)==0)
%                 figure(1)
%                 clf
%                 plot(xx,sum(Mnew*dalpha),'-')
%                 hold on
%                 plot(xx,sum(vnew*dalpha),'-')
%                 drawnow
%             end
        end
        V(:,:,iy,ia) = vnew;
    end
    save('Evol_lineages','V','Mnew')
end

% %% 
% load('Evol_lineages.mat')
% %%
% M_x = sum(Mnew*dalpha);
% V_x = permute(sum(V*dalpha,1),[3,4,2,1]);
% V_x = reshape(V_x,Na*Ny,Nx);
% Vx = cumsum(V_x,1);
% figure(1)
% clf
% plot(xx,M_x,'--')
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
% %% Figure Proportion
% Prop = V_x./M_x;
% prop= cumsum(Prop);
% figure(2)
% clf
% plot(xx,prop)



%% Figures
% Color = get(gca,'colororder');
% 
% figure(11)
% clf
% % yyaxis left
% scatter(BETA,C_simu(1,:),'filled')
% hold on
% % plot(BETA,c1_star,'-.')
% % plot(BETA,c2_star,'--')
% % plot(BETA,c3_star,':')
% % plot(BETA,c4_star,'d')
% plot(BETA,c5_star,'--','Color',Color(1,:))
% % plot(BETA,c6_star,'*')
% ylabel('Spreading speed (AMF and plant)','Interpreter','latex','FontSize',20)
% hold on
% % scatter(BETA,C_simu(2,:),'filled')
% % line([BETA(1),BETA(end)],[0,0])
% ylim([0,0.15])
% % yyaxis right
% % hold on
% % plot(BETA,sqrt(PSTAR))
% % % plot(BETA,MSTAR)
% % ylabel('Plant and AMF biomass','Interpreter','latex')
% xlabel('Carbon suupply rate ($\beta$)','Interpreter','latex','FontSize',20)
% 
% figure(2)
% clf
% yyaxis left
% scatter(BETA,ABARSTAR,'filled')
% hold on
% plot(BETA,ABARSTAR_approx)
% yyaxis right
% hold on
% scatter(BETA,PSTAR,'d','filled')
% plot(BETA,PSTAR_approx)
% scatter(BETA,MSTAR,'filled')
% plot(BETA,MSTAR_approx)
% 
% 
% 
% 
