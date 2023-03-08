clear all
close all
% format LongE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Michorizer model with frequency dependence   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time parameter
Tf = 70;
%% Parameter of the model
%% Space
xmin = -10;
xmax = 100;
dx = 0.05;
xx = xmin:dx:xmax;
Nx = length(xx);

% Diffusion matrix
e = ones(Nx,1);
I  = spdiags(e,0,Nx,Nx);
Ad = spdiags([e -2*e e],-1:1,Nx,Nx);
Ad(1,1) = -1;
Ad(end,end) = -1;
Ad = Ad/(dx^2);
dp = 0.1;
dm = 0.1;

% Advection matrix
Ac =  spdiags([-e e],0:1,Nx,Nx);
Ac(end,end) = 0;

global q d mui mup
q = 6;

mup = 0.3; % 1/100
mui = 0.3; % 1/20

d = 1.2;
%% Growth rate plant
RP = [0,0.2,0.8];
Nrp = length(RP);
%% Trait alpha AMF
ALPHA = [0.4,0.7,0.9];
Nalpha= length(ALPHA);
%% Trait Beta plant
Nbeta = 20;
BETA = zeros(Nalpha,Nbeta);

C_simu = zeros(Nalpha,Nbeta);
C_approx = zeros(Nalpha,Nbeta);

for irp = 1:Nrp
    rp = RP(irp);
    
    for ia = 1:Nalpha
        %% Alpha trait
        alpha = ALPHA(ia);
        
        %% Beta trait
        Pstar  =@(beta) (q+1+sqrt( (q-1)^2 -4*q*mup*mui./(beta.^2)) )./(2*(beta+mup*mui./beta));
        betamax = fzero(@(x) alpha*Pstar(x)-d,13);
        % betamin = fzero(@(x) alpha*Pstar(x)-d,0.1);
        betamin = sqrt(4*q*mup*mui./((q-1).^2));
        % betamax = alpha/d;
        %max(2.5,betamin);
        BETA(ia,:) = linspace(betamin+0.05,betamax-0.1,Nbeta);
        parfor ib = 1:Nbeta
            beta = BETA(ia,ib);
            % Boundary condition
            %% rp =0
            %         p_star = alpha*Pstar(beta)-d;
            %         m_star = p_star*(beta-1./Pstar(beta))/mui;
            %% rp >0
            RP = [-alpha*(beta^2+mui*mup), rp*mui+beta*alpha*(q+1)+d*(beta^2+mui*mup),...
                -beta*(q+1)*d - alpha*q , q*d];
            X = alpha*roots(RP)-d;
            p_star = max(X(X>0));
            m_star = p_star*(beta-alpha./(p_star+d))/mui;
            
            %% 1 equal speed
            % gam = sqrt((beta-alpha/(d))./(q*alpha/(d)-beta));
            % rd = sqrt((q*alpha/(d)-beta).*(beta-alpha/(d))) -mup;
            %% 2 proportion m/p
            u_star = (beta-alpha/d + mup)./(q*alpha/d-beta+mui);
            rd = (q*alpha/d-beta)*u_star-mup;
            c_star = sqrt(dp*rd/2);
            %% 3rd approx
            G = [q*alpha/d-beta, mui*m_star./p_star+mup,-m_star./p_star.*(beta-alpha/d)];
            GAM = max(roots(G));
            rd_5 = (q*alpha/d-beta).*GAM-mup;
            c5_star = sqrt(dp.*p_star.*real(rd_5)/2);
            C_approx(ia,ib,irp) = c5_star;
            % rd_minus = (Q*alpha/d-beta).*(beta-alpha/d) -mup;
            % c_minus_star = sqrt(D_p*rd_minus/2);
            %         else
            %             p_star = 0;
            %             m_star = 0;
            %             C_approx(ia,ib) = 0;
            %         end
            fp =@(p,m) p.*(rp+(q*alpha./(p+d) -beta).*m - mup*p);
            fm =@(p,m) m.*((beta-alpha./(p+d) ).*p - mui*m);
            
            
            
            %% Solver PDE
            it = 1;
            dt = 0.05;
            tt = 0:dt:Tf;
            Nt = length(tt);
            
            PP = zeros(Nt,Nx); MM = zeros(Nt,Nx);
            P0 =  p_star.*(xx<0);
            M0 = m_star.*(xx<0);
            
            Pnew   = P0';  PP(1,:) = P0;
            Mnew  = M0';   MM(1,:) = M0;
            
            Xt = zeros(1,Nt);
            m_speed = m_star/2;
            Ic =  max(1,sum(Mnew<m_speed)); xt = xx(Ic); Xt(it) = xt;
            err=1;
            
            % figure(1)
            % clf
            % hold on
            % plot(xx,Pnew)
            % plot(xx,Mnew)
            % drawnow
            
            while (it<Nt)%&&(err>1e-4)
                Pold = Pnew;Mold = Mnew;
                %% Ref frame
                Pnew   = (I-dt*dp*Ad)\(Pold + dt*fp(Pold,Mold));
                Mnew  = (I-dt*dm*Ad)\(Mold + dt*fm(Pold,Mold));
                
                %% Moving frame at speed c_cstar
                %     Pnew   = (I-dt*dp*Ad  - dt*c_star*Ac/dx)\(Pold + dt*fp(Pold,Mold)  );% + dt*BP/(dx^2) +dt*Bp/dx);
                %     Mnew  = (I-dt*dm*Ad - dt*c_star*Ac/dx)\(Mold + dt*fm(Pold,Mold) );% + dt*BMc/(dx^2));
                Ic =  (Mnew<m_speed);
                xt = min([xx(Ic),xx(end)]);
                it = it + 1;
                Xt(it)   = xt;
                PP(it,:) = Pnew';
                MM(it,:) = Mnew';
                %     err = max(abs(Pold-Pnew));
                %% Figure
                % if (mod(it,10)==0)
                % clf
                % hold on
                % plot(xx,Pnew,'--')
                % plot(xx,Mnew)
                % drawnow
                % end
            end
            tt = 0:dt:Tf;
            c_simu = mean(Xt(end-50:end)/(tt(end-50:end)));
            C_simu(ia,ib,irp) = c_simu;
        end
    end
end
%% figure speed
Color = get(gca,'colororder');
marker = {'o','d','^','v'};
lines = {'-','--','-.'};
cstar_p = 2*sqrt(dp*RP);

figure(2)
clf
hold on
p = zeros(1,Nrp);
for irp = 1:Nrp

for ia = 1:Nalpha
    p(irp) = scatter(BETA(ia,:),C_simu(ia,:,irp),'filled','MarkerFaceColor',Color(ia,:),...
            'Marker',marker{irp},'DisplayName', ['$r_p= $', num2str(RP(irp))]);
    hold on
     plot(BETA(ia,:),C_approx(ia,:,irp),'LineStyle',lines{irp},'color',Color(ia,:))
end
line([0,5],[cstar_p(irp),cstar_p(irp)],'color','k','linestyle',lines{irp})
end
xlim([betamin,betamax])
hl = legend(p);
set(hl, 'Interpreter','latex', 'FontSize',16)
xlabel('Carbon suupply rate ($\beta$)','Interpreter','latex','fontsize',20)
ylabel('Spreading speed (AMF and plant)','Interpreter','latex','fontsize',20)

%% Relative spreading speed
% figure(2)
% clf
% hold on
% p = zeros(1,Nrp);
% for irp = 1:Nrp
% DC_simu = C_simu(:,:,irp) - cstar_p(irp);
% DC_approx = C_approx(:,:,irp) - cstar_p(irp);
% for ia = 1:Nalpha
%     p(irp) = scatter(BETA(ia,:),DC_simu(ia,:),'filled','MarkerFaceColor',Color(ia,:),...
%             'Marker',marker{irp},'DisplayName', ['$r_p= $', num2str(RP(irp))]);
% %     hold on
% %      plot(BETA(ia,:),DC_approx(ia,:),'LineStyle',lines{irp},'color',Color(ia,:))
% end
% end
% line([0,5],[0,0],'color','k')
% xlim([betamin,betamax])
% hl = legend(p);
% set(hl, 'Interpreter','latex', 'FontSize',16)
% xlabel('Carbon suupply rate ($\beta$)','Interpreter','latex','fontsize',20)
% ylabel('Spreading speed (AMF and plant)','Interpreter','latex','fontsize',20)

% figure(2)
% hold on
% plot(tt,Xt,'--')
% plot(tt,c_star*tt)
% plot(tt,c5_star*tt)



