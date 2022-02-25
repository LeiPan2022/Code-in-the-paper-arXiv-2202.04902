%This program plots the average Sz for both multi-magnon states and eigenstates in the PPXPP model and the overlap between them.
%Author: Lei Pan, Email:panlei@mail.tsinghua.edu.cn

L=18;

[Ener_k,VV_k,Sz_average_Mom]=PPXPP_Ener_Sz_Momentum(L);
[Max_Overlap2,Max_index2,Ener_Multi_Magnon_PPXPP,Sz_average_Magnons]=Multi_Magnon_Ener_PPXPP(L);


figure
subplot(2,1,1)
hold on; plot(1:15,Max_Overlap2,'.','LineWidth',3,'MarkerSize',60)
xlabel('Magnon number n','fontsize',40);
ylabel('Overlap','fontsize',40);
set(gca,'linewidth',3);
set(gca,'FontSize',36);
box on;    
set(gca, 'FontName', 'Times New Roman'); 

subplot(2,1,2)
hold on;
plot(Ener_k,Sz_average_Mom,'.','LineWidth',3,'MarkerSize',60)

  hold on;
  plot(Ener_k(Max_index2(1:15)),Sz_average_Mom(Max_index2(1:15)),'o','LineWidth',3,'MarkerSize',23)
 hold on;
 plot(Ener_Multi_Magnon_PPXPP,Sz_average_Magnons,'X','LineWidth',3,'MarkerSize',23)
 xlabel('Energy','fontsize',40);
ylabel('S_z','fontsize',40);
set(gca,'linewidth',3);
set(gca,'FontSize',36);
box on;  
set(gca, 'FontName', 'Times New Roman'); 

