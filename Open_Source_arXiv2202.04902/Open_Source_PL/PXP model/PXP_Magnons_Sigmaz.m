function [Sigmaz_Magnons_PXP]=PXP_Magnons_Sigmaz(L,Pi_Magnons_vec)


 [Dim_sub,vector,E_PXP,VV,find_index,index,~]=PXP_Ham(L);
% [Overlap_PXP,Overlap_PXP_index,Pi_Magnons_vec]=PXP_Overlap_Magnons(L)

Sz_Magnons_PXP=zeros(L,L/2+1);
for jj=1:L/2+1
    
    for kk=1:L
       
        Sz_Magnons_PXP(kk,jj)=(2*vector(:,kk)-1)'*(Pi_Magnons_vec(:,jj).*conj(Pi_Magnons_vec(:,jj)));       
               
    end
    
end

Sigmaz_Magnons_PXP(1,:)=sum(Sz_Magnons_PXP)/L;

end