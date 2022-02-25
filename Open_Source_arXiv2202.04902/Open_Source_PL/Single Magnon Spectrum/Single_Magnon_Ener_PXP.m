function [Ener_single_Magnon_PXP]=Single_Magnon_Ener_PXP(L)

uu=1;  
Dim=2^L;
[Dim_sub,vector,E_PXP,VV,find_index,index,Ham_PXP]=PXP_Ham(L);
vv=Find_v_PXP(L);
phi_uv=cell(1,L+1);

phi_uv{1}=1;
phi_unit=[uu;vv]; phi_flip=[-vv;uu];

BCS_Basis=zeros(Dim_sub,1+L); BCS_Basis0=zeros(Dim,1+L);


for ll=1:L
    
    phi_uv{ll+1}=kron(phi_uv{ll},phi_unit);
    
end
    
BCS_Basis0(:,1)=phi_uv{L+1};

BCS_Basis(:,1)=BCS_Basis0(index,1);

BCS_Basis(:,1)=BCS_Basis(:,1)/norm(BCS_Basis(:,1)); 



for kk=1:L

nk=(kk-1);
Momentum=exp(1i*2*pi*nk/L);
for ll=1:L
    sign=(Momentum)^(ll);
    
    BCS_Basis0(:,kk+1)=BCS_Basis0(:,kk+1)+sign*kron(kron(phi_uv{ll},phi_flip),phi_uv{L+1-ll});
   % BCS_Basis01(:,kk+1)=BCS_Basis01(:,kk+1)+sign*BCS_Basis3(:,ll+1);
end

BCS_Basis(:,kk+1)=BCS_Basis0(index,kk+1);
BCS_Basis(:,kk+1)=BCS_Basis(:,kk+1)/norm(BCS_Basis(:,kk+1)); 

end

BCS_Basis2=zeros(Dim_sub,L+1);
BCS_Basis2(:,1)=BCS_Basis(:,1);
for ii=2:L+1
for jj=1:ii-1
BCS_Basis2(:,ii)=BCS_Basis2(:,ii)-(BCS_Basis(:,ii)'*BCS_Basis2(:,jj))/(BCS_Basis2(:,jj)'*BCS_Basis2(:,jj))*BCS_Basis2(:,jj);
end
BCS_Basis2(:,ii)=BCS_Basis2(:,ii)+BCS_Basis(:,ii);
end

for k=1:L+1
BCS_Basis2(:,k)=BCS_Basis2(:,k)/norm(BCS_Basis2(:,k));
end

%Single-magnon
for mm=1:L
    mm
    for nn=1:L
    
        %Ham_magnon(mm,nn)=BCS_Basis(:,mm+1)'*Ham_PXP*BCS_Basis(:,nn+1);
        Ham_magnon(mm,nn)=BCS_Basis2(:,mm+1)'*Ham_PXP*BCS_Basis2(:,nn+1);
    end
end


[V_m,Ener_m]=eig(Ham_magnon);
[DD_m,J_m]=sort(real(diag(Ener_m)),'ascend');
  
Ener_single_Magnon_PXP=[diag(real(Ham_magnon));Ham_magnon(1,1)]-E_PXP(1);
figure
hold on;
plot(([0:1:L]/L),Ener_single_Magnon_PXP,'.','LineWidth',3,'MarkerSize',60)%,'MarkerEdgeColor','m')
xlabel('$k/\pi$','fontsize',40);
ylabel('Energy','fontsize',40);
set(gca,'linewidth',3);
set(gca,'FontSize',36);
box on;    
set(gca, 'FontName', 'Times New Roman'); 

end