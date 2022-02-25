%This program plots the single-magnon excitation spectrum in the PXP model and the PPXPP model.
%Author: Lei Pan, Email:panlei@mail.tsinghua.edu.cn

 

L=18;
[Dim_sub_PPXPP,vector_PPXPP,E_PPXPP,VV_PPXPP,find_index_PPXPP,index_PPXPP,Ham_PPXPP]=PPXPP_Ham(L);
[Ener_single_Magnon_PXP]=Single_Magnon_Ener_PXP(L);
[Ener_single_Magnon_PPXPP]=Single_Magnon_Ener_PPXPP(L);
 
figure
 
plot(([0:1:L]/L),Ener_single_Magnon_PXP,'.','LineWidth',3,'MarkerSize',60) 
hold on;
plot(([0:1:L]/L),Ener_single_Magnon_PPXPP,'o','LineWidth',3,'MarkerSize',25) 
xlabel('$k/2\pi$','fontsize',40);
ylabel('Energy','fontsize',40);
set(gca,'linewidth',3);
set(gca,'FontSize',36);
box on;    
set(gca, 'FontName', 'Times New Roman'); 

xi=0:0.01:1;
y1=spline(([0:1:L]/L),Ener_single_Magnon_PXP,xi);
y2=spline(([0:1:L]/L),Ener_single_Magnon_PPXPP,xi);
plot(xi,y1,'-','LineWidth',4,'Color',[0 0.447 0.741]); 
hold on
plot(xi,y2,'--','LineWidth',4,'Color',[0.850980401039124 0.325490206480026 0.0980392172932625]);
