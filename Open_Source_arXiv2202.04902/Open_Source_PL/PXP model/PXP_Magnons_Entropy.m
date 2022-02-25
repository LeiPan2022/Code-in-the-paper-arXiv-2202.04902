function [Entropy_Magnons_PXP]=PXP_Magnons_Entropy(L,Pi_Magnons_vec)


% [Dim_sub,vector,E_PXP,VV,find_index,index,~]=PXP_Ham(L);
% [Overlap_PXP,Overlap_PXP_index,Pi_Magnons_vec]=PXP_Overlap_Magnons(L)
[Dim_sub,vector,Dim_sub2,vector2]=Vector_PXP(L);

num_site=L/2; 
Entropy_Magnons_PXP=zeros(L/2+1,1);
 
for kk=1:L/2+1
    
    kk
    rho0=zeros(Dim_sub2,Dim_sub2); 
  state=Pi_Magnons_vec(:,kk);

for ii=1:Dim_sub2
    for jj=ii:Dim_sub2      
        
        find_index_ii=find(sum(abs(vector2(ii,:)-vector(:,1:num_site)),2)==0);  
        find_index_jj=find(sum(abs(vector2(jj,:)-vector(:,1:num_site)),2)==0);
        
        [C,index_ii,index_jj]=intersect(vector(find_index_ii,num_site+1:end),vector(find_index_jj,num_site+1:end),'rows');
        
        
        
        rho0(ii,jj)=state(find_index_jj(index_jj),1)'*state(find_index_ii(index_ii),1); 
        rho0(jj,ii)=rho0(ii,jj)';
    end
end
Rho_temp=rho0;


Entropy_Magnons_PXP(kk,1)=-trace(Rho_temp*logm(Rho_temp+0.000000000001));
warning('off')
end
 

Entropy_Magnons_PXP=real(Entropy_Magnons_PXP);

end