function [Sigmaz_Magnons_PXP2]=PXP_Magnons_Sigmaz2(L,Pi_Magnons_vec2)


 [Dim_sub,vector,E_PXP,VV,find_index,index,~]=PXP_Ham(L);
% [Overlap_PXP,Overlap_PXP_index,Pi_Magnons_vec]=PXP_Overlap_Magnons(L)

Sz_Magnons_PXP2=zeros(L,L/2+1);
for jj=1:L/2+1
    
    for kk=1:L
       
        Sz_Magnons_PXP2(kk,jj)=(2*vector(:,kk)-1)'*(Pi_Magnons_vec2(:,jj).*conj(Pi_Magnons_vec2(:,jj)));       
               
    end
    
end

Sigmaz_Magnons_PXP2(1,:)=sum(Sz_Magnons_PXP2)/L;
%Sigmaz_Magnons_PXP2=fliplr(Sigmaz_Magnons_PXP2);


end



