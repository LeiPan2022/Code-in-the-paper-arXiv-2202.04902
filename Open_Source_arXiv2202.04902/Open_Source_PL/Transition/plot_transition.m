%This program plots the optimal v for magnon wave function, the overlap between magnon wave function and the ground state obtained by ED, the gap of magnon excitation and the domain wall density and variance in the PXP model with finite detuning.
%Author: Lei Pan, Email:panlei@mail.tsinghua.edu.cn


 L=12; 
 Detuning=-0.5:0.05:1.5;

for bb=1:length(Detuning)
    bb
    B=Detuning(bb); 

    BB(bb)=B;
[vv_B,Gap_B,Overlap_B,Ener_B,Groundstate_B]=PXP_Finite_detuning(L,B);


Gap_B2(bb)=Gap_B;
vv_B2(bb)=vv_B;
Overlap_B2(bb)=Overlap_B;

[Domain_wall_mean,Domain_wall_flu]=PXP_Domain_Wall_Density(L,Groundstate_B);

Domain_wall_mean_B(bb)=Domain_wall_mean;
Domain_wall_flu_B(bb)=Domain_wall_flu;

end


 
figure
hold on;
subplot(2,2,1)

plot(BB,vv_B2,'.-','LineWidth',3,'MarkerSize',60) 
hold on; 
 
xlabel('$\Delta$','Interpreter','latex','fontsize',40);
ylabel('v','fontsize',40);
set(gca,'linewidth',3);
set(gca,'FontSize',36);
set(gca, 'FontName', 'Times New Roman'); 
box on;  


subplot(2,2,2)
hold on; 
plot(BB,abs(Overlap_B2),'.-','LineWidth',3,'MarkerSize',60) 

xlabel('$\Delta$','Interpreter','latex','fontsize',40);
ylabel('Overlap','fontsize',40);set(gca,'linewidth',3);
set(gca,'FontSize',36);
set(gca, 'FontName', 'Times New Roman'); 
box on;  


subplot(2,2,3)
hold on;
plot(BB,Gap_B2,'.-','LineWidth',3,'MarkerSize',60) 
 
xlabel('$\Delta$','Interpreter','latex','fontsize',40);
ylabel('Gap','fontsize',40);set(gca,'linewidth',3);
set(gca,'FontSize',36);
set(gca, 'FontName', 'Times New Roman'); 
box on;  



subplot(2,2,4)
[AX,H1,H2] = plotyy(BB,Domain_wall_mean_B,BB,Domain_wall_flu_B,'plot');
set(AX,'FontSize',40)
set(AX,'linewidth',4)

% set(AX(1),'XColor','k','YColor','b');
% set(AX(2),'XColor','k','YColor','r');

set(AX(1),'XColor','k','YColor','[0    0.4470    0.7410]');
set(AX(2),'XColor','k','YColor','[0.8500    0.3250    0.0980]');

HH1=get(AX(1),'Ylabel');
set(HH1,'String','Domain-wall density');
set(HH1,'color','[0    0.4470    0.7410]');
set(gca,'linewidth',3); 


HH2=get(AX(2),'Ylabel');
set(HH2,'String','Variance');
set(HH2,'color','[0.8500    0.3250    0.0980]');
set(gca,'linewidth',4); 

set(H1,'LineStyle','-','LineWidth',6);
%set(H1,'color','b');
set(H2,'LineStyle','--','LineWidth',6);
%set(H2,'color','r');
rgb_value1=get(H1,'Color');
rgb_value2=get(H2,'Color');
%set(AX,'FontSize',32)
legend([H1,H2],{'Mean';'Flutuation'});
xlabel('$\Delta$','Interpreter','latex','fontsize',40);
%title('Labeling plotyy');
% set(gca,'linewidth',3);
% set(gca,'FontSize',32);
box on;
set(gca, 'FontName', 'Times New Roman');
%axis square
set(gca, 'XTick', [-0.5 0 0.5 1 1.5])
set(gca,'XLim',[-0.5 1.5]); 
set(gca,'TickLabelInterpreter', 'latex');

 































