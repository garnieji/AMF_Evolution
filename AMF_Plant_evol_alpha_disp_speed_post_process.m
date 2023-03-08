clear all
close all
% format LongE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Michorizer model with frequency dependence and evolution of alpha no space   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time parameter
Color = get(gca,'colororder');
Data = {'data_speed_low_plant.mat','data_speed.mat','data_speed_low_AMF.mat'};
Nd = length(Data);
p =zeros(1,Nd);
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
    p(idata) = scatter(BETA,C_simu(1,:),'filled','MarkerFaceColor',Color(idata,:),'DisplayName',['$\frac{D_p}{D_m}= $',num2str(DD)]);
    hold on
    plot(BETA,c5_star,'-','Color',Color(idata,:))
end

ylabel('Spreading speed (AMF and plant)','Interpreter','latex','FontSize',20)
ylim([0,0.15])
hl = legend(p);
set(hl, 'Interpreter','latex', 'FontSize',16)
xlabel('Carbon suupply rate ($\beta$)','Interpreter','latex','FontSize',20)

figure(2)
yyaxis left
scatter(BETA,ABARSTAR,'filled','MarkerFaceColor',Color(1,:))
hold on
plot(BETA,ABARSTAR_approx,'Color',Color(1,:))
yyaxis right
hold on
scatter(BETA,PSTAR,'d','filled')
plot(BETA,PSTAR_approx)
scatter(BETA,MSTAR,'filled')
plot(BETA,MSTAR_approx)


