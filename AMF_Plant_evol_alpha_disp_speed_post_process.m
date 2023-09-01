clear all
close all
% format LongE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Michorizer model with frequency dependence and evolution of alpha no space   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time parameter
Color = get(gca,'colororder');
% Data = {'data_speed_1_low_plant.mat','data_speed_0.mat','data_speed_1_low_AMF.mat'};
DATA = {{'data_speed_slow_plant.mat','data_speed_same_rp_0.mat','data_speed_slow_AMF.mat'};...
        {'data_speed_same_rp_0.mat','data_speed_same_rp_01.mat','data_speed_same_rp_02.mat','data_speed_same_rp_08.mat'}};

%% Figure Diff diffusivity
Data = DATA{1};
Nd = length(Data);
p = zeros(1,Nd);
figure(1)
clf
for idata = 1:Nd
    load(Data{idata})
    DD = D_p/D_m;
    
    %% Speed approx
    G = [Q*ABARSTAR_approx/d-BETA;...
        mui*D_m/D_p*MSTAR_approx./PSTAR_approx+mup;...
        -D_m/D_p*MSTAR_approx./PSTAR_approx.*(BETA-ABARSTAR_approx/d)]';
    GAM = zeros(1,Nbeta);
    for ib = 1:Nbeta
        GAM(ib) = max(roots(G(ib,:)));
    end
    rd_5 = (Q*ABARSTAR_approx/d-BETA).*GAM-mup;
    c5_star = sqrt(D_p.*PSTAR_approx.*rd_5/2);
    
    %% Figures  
    figure(1)
    % yyaxis left
    if (idata==1)
    semilogx(BETA,c5_star,'-','Color',Color(idata,:))
    end
    p(idata) = scatter(BETA,C_simu(1,:),'filled','MarkerFaceColor',Color(idata,:),'DisplayName',['$\frac{D_p}{D_m}= $',num2str(DD)]);
    hold on
    semilogx(BETA,C_simu_approx(1,:),'Color',Color(idata,:))
end
set(gca,'Xscale','log','XTick',[0.4,1,4,10,12],'XTickLabel',{'0.4','$10^0$','4','$10^1$','12'},'TickLabelInterpreter','latex','FontSize',14)
ylabel('Speed of spread','Interpreter','latex','FontSize',20)
axis([betamin,betamax,0,0.18])
hl = legend(p);
set(hl, 'Interpreter','latex', 'FontSize',16)
xlabel('Carbon supply rate,  $\beta$ (log scale)','Interpreter','latex','FontSize',20)

%% Figure Diff rp 
Data = DATA{2};
Nd = length(Data);
p = zeros(1,Nd);

%% Speed
figure(2)
clf
for idata = 1:Nd
    load(Data{idata})
    Cp = 2*sqrt(D_p*rp);
    DD = D_p/D_m;
    
    %% Speed approx
    G = [Q*ABARSTAR_approx/d-BETA;...
        mui*D_m/D_p*MSTAR_approx./PSTAR_approx+mup;...
        -D_m/D_p*MSTAR_approx./PSTAR_approx.*(BETA-ABARSTAR_approx/d)]';
    GAM = zeros(1,Nbeta);
    for ib = 1:Nbeta
        GAM(ib) = max(roots(G(ib,:)));
    end
    rd_5 = (Q*ABARSTAR_approx/d-BETA).*GAM-mup;
    c5_star = sqrt(D_p.*PSTAR_approx.*rd_5/2);
    
    Diff_C = C_simu(1,:)-Cp;  
    Diff_C_approx = C_simu_approx(1,:)-Cp;  
    %% Figures  
    figure(2)
    % yyaxis left
    if (idata==1)
    semilogx(BETA,c5_star,'-','Color',Color(idata,:))
    end
    p(idata) = scatter(BETA,C_simu(1,:),'filled','MarkerFaceColor',Color(idata,:),'DisplayName',['$r_p= $',num2str(rp)]);
    hold on
    semilogx(BETA,C_simu_approx(1,:),'Color',Color(idata,:))
    line([betamin,betamax],[Cp,Cp],'color',Color(idata,:),'linewidth',2,'linestyle','--')
end
set(gca,'Xscale','log','XTick',[0.4,1,4,10,12],'XTickLabel',{'0.4','$10^0$','4','$10^1$','12'},'TickLabelInterpreter','latex','FontSize',14)
axis([betamin,betamax,0,2*sqrt(D_p*rp)+0.05])
ylabel('Speed of spread','Interpreter','latex','FontSize',20)
hl = legend(p);
set(hl, 'Interpreter','latex', 'FontSize',16)
xlabel('Carbon supply rate, $\beta$ (log scale)','Interpreter','latex','FontSize',20)


%% Diff speed
figure(21)
clf
for idata = 1:Nd
    load(Data{idata})
    Cp = 2*sqrt(D_p*rp);
    DD = D_p/D_m;
    
    %% Speed approx
    G = [Q*ABARSTAR_approx/d-BETA;...
        mui*D_m/D_p*MSTAR_approx./PSTAR_approx+mup;...
        -D_m/D_p*MSTAR_approx./PSTAR_approx.*(BETA-ABARSTAR_approx/d)]';
    GAM = zeros(1,Nbeta);
    for ib = 1:Nbeta
        GAM(ib) = max(roots(G(ib,:)));
    end
    rd_5 = (Q*ABARSTAR_approx/d-BETA).*GAM-mup;
    c5_star = sqrt(D_p.*PSTAR_approx.*rd_5/2);
    
    Diff_C = C_simu(1,:)-Cp;  
    Diff_C_approx = C_simu_approx(1,:)-Cp;  
    %% Figures  
    figure(21)
    % yyaxis left
    if (idata==1)
    semilogx(BETA,c5_star,'-','Color',Color(idata,:))
    end
%     p(idata) = scatter(BETA,C_simu(1,:),'filled','MarkerFaceColor',Color(idata,:),'DisplayName',['$r_p= $',num2str(rp)]);
    p(idata) = scatter(BETA,Diff_C,'filled','MarkerFaceColor',Color(idata,:),'DisplayName',['$r_p= $',num2str(rp)]);

    hold on
%     semilogx(BETA,C_simu_approx(1,:),'Color',Color(idata,:))
    semilogx(BETA,Diff_C_approx,'Color',Color(idata,:))

end
line([betamin,betamax],[0,0],'color','k')
set(gca,'Xscale','log','XTick',[0.4,1,4,10,12],'XTickLabel',{'0.4','$10^0$','4','$10^1$','12'},'TickLabelInterpreter','latex','FontSize',14)
ylabel({'Difference of spreading speed', 'with and without AMF community'},'Interpreter','latex','FontSize',20)
% axis([betamin,betamax,0,0.18])
hl = legend(p);
set(hl, 'Interpreter','latex', 'FontSize',16)
xlabel('Carbon supply rate, $\beta$ (log scale)','Interpreter','latex','FontSize',20)


%% Optimal BETA

load('data_speed_same_rp_0.mat')

[P_max,Ipmax] = max(PSTAR);
Beta_P_max = BETA(Ipmax);
[M_max,Immax] = max(MSTAR);
Beta_M_max = BETA(Immax);
[C_max,Icmax] = max(C_simu(1,:));
Beta_C_max = BETA(Icmax);
figure(3)
clf
yyaxis left
semilogx(BETA,C_simu_approx(1,:))
hold on
scatter(BETA,C_simu(1,:),'filled','MarkerFaceColor',Color(1,:))
% scatter(BETA,ABARSTAR,'filled','MarkerFaceColor',Color(1,:))
% hold on
% plot(BETA,ABARSTAR_approx,'Color',Color(1,:))

yyaxis right
hold on
scatter(BETA,PSTAR,'d','filled')
plot(BETA,PSTAR_approx)
scatter(BETA,MSTAR,'filled')
plot(BETA,MSTAR_approx)


% 
% 
