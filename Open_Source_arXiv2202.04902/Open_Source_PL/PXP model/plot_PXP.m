%This program plots the overlap between multi-magnon states and scar states, and also the average Sz and bipartite entanglement entropy for both multi-magnon states and scar states in the PXP model 
%Author: Lei Pan, Email:panlei@mail.tsinghua.edu.cn

 

L=12;

[Dim_sub,vector,E_PXP,VV,find_index,index,Ham_PXP]=PXP_Ham(L);
[Dim_sub,vector,E_PXP,Sz_PXP,find_index,index]=PXP_Ener_Sz(L);
[Dim_sub,vector,E_PXP,Entropy_PXP,find_index,index]=PXP_Ener_Entropy(L);
[Overlap_PXP,Overlap_PXP_index,Pi_Magnons_vec,Ener_Magnons]=PXP_Overlap_Magnons(L);
[Overlap_PXP2,Overlap_PXP_index2,Pi_Magnons_vec2,Ener_Magnons2]=PXP_Overlap_Magnons2(L);
[Sigmaz_Magnons_PXP]=PXP_Magnons_Sigmaz(L,Pi_Magnons_vec);
[Sigmaz_Magnons_PXP2]=PXP_Magnons_Sigmaz(L,Pi_Magnons_vec2);

[Entropy_Magnons_PXP]=PXP_Magnons_Entropy(L,Pi_Magnons_vec);
[Entropy_Magnons_PXP2]=PXP_Magnons_Entropy2(L,Pi_Magnons_vec2);




figure
hold on;
subplot(3,1,1)

plot(1:L/2+1,Overlap_PXP,'.','LineWidth',3,'MarkerSize',60) 
hold on; 
plot((L+1):-1:L/2+1,Overlap_PXP2,'.','LineWidth',3,'MarkerSize',60) 
xlabel('Scar index n','fontsize',40);
ylabel('Overlap','fontsize',40);
set(gca,'linewidth',3);
set(gca,'FontSize',36);
set(gca, 'FontName', 'Times New Roman'); 
box on;  


subplot(3,1,2)

plot(E_PXP,Sz_PXP/2,'.','LineWidth',3,'MarkerSize',60) 
hold on; 
plot(Ener_Magnons,Sigmaz_Magnons_PXP/2,'X','LineWidth',3,'MarkerSize',25) 
hold on; 
plot(Ener_Magnons2,Sigmaz_Magnons_PXP2/2,'X','LineWidth',3,'MarkerSize',25) 
xlabel('Energy','fontsize',40);
ylabel('$-\langle S_z\rangle$','Interpreter','latex','fontsize',40);
set(gca,'linewidth',3);
set(gca,'FontSize',36);
set(gca, 'FontName', 'Times New Roman'); 
box on;  


subplot(3,1,3)

plot(E_PXP,Entropy_PXP,'.','LineWidth',3,'MarkerSize',60) 
hold on;
plot(Ener_Magnons,Entropy_Magnons_PXP,'X','LineWidth',3,'MarkerSize',25) 
hold on; 
plot(Ener_Magnons2,Entropy_Magnons_PXP2,'X','LineWidth',3,'MarkerSize',25) 
xlabel('Energy','fontsize',40);
ylabel('$S_v$','Interpreter','latex','fontsize',40);
set(gca,'linewidth',3);
set(gca,'FontSize',36);
set(gca, 'FontName', 'Times New Roman'); 
box on;  

