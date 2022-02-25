function [Dim_sub,vector,E_PXP,Sz_PXP,find_index,index]=PXP_Ener_Sz(L)


[Dim_sub,vector,E_PXP,VV,find_index,index,~]=PXP_Ham(L);


Sz_1average=zeros(L,Dim_sub);
for jj=1:Dim_sub
    
    for kk=1:L
       
        Sz_1average(kk,jj)=(2*vector(:,kk)-1)'*(VV(:,jj).*conj(VV(:,jj)));       
               
    end
    
end

Sz_PXP(1,:)=sum(Sz_1average)/L;

end