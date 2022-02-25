function [Dim_sub_PPXPP,vector_PPXPP,E_PPXPP,Sz_PPXPP,find_index_PPXPP,index_PPXPP]=PPXPP_Ener_Sz(L)


[Dim_sub_PPXPP,vector_PPXPP,E_PPXPP,VV_PPXPP,find_index_PPXPP,index_PPXPP,Ham_PPXPP]=PPXPP_Ham(L);


Sz_1average=zeros(L,Dim_sub_PPXPP);
for jj=1:Dim_sub
    
    for kk=1:L
       
        Sz_1average(kk,jj)=(2*vector_PPXPP(:,kk)-1)'*(VV_PPXPP(:,jj).*conj(VV_PPXPP(:,jj)));       
               
    end
    
end

Sz_PXP(1,:)=sum(Sz_1average)/L;

end