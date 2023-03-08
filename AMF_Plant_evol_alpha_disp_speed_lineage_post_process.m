clear all
close all
% format LongE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Michorizer model with frequency dependence                           %
%        evolution of alpha spreading in space                            %
%          LINEAGES                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
load('Evol_lineages.mat')
%%
M_x = sum(Mnew*dalpha);
V_x = permute(sum(V*dalpha,1),[3,4,2,1]);
V_x = reshape(V_x,Na*Ny,Nx);
Vx = cumsum(V_x,1);
figure(1)
clf
plot(xx,M_x,'--')
hold on
% plot(xx,sum(Minit*dalpha))
% if (iy<Ny)
% plot(xx,sum(Minit.*((xx<yy(iy+1)).*(xx>=yy(iy)))*dalpha))
% else
% plot(xx,sum(Minit.*((xx>yy(iy)))*dalpha))
% end 
for ii = 1:Ny*Na
        plot(xx,Vx(ii,:),'-')
end
%% Figure Proportion
Prop = V_x./M_x;
prop= cumsum(Prop);
figure(2)
clf
plot(xx,prop)