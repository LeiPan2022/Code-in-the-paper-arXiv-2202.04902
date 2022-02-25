function [Dim_sub,vector,E_PXP,Entropy_PXP,find_index,index]=PXP_Ener_Entropy(L)


[Dim_sub,vector,E_PXP,VV,find_index,index,~]=PXP_Ham(L);


num_site=L/2;

 
temp2=zeros(1,num_site); temp3=zeros(1,num_site); vec2=zeros(1,num_site);
tic; cc=0;
for jj=1:2^num_site
    
obj=dec2base(jj-1,2,num_site);  
temp2=(obj=='1')*1+(obj=='0')*0;
temp3(2:end)=temp2(1:end-1); 
 
 
kk2=sum(temp2.*temp3);  
   if kk2==0
       cc=cc+1;
      vec2(cc,:)=temp2(1,:); 
      
   end
end

toc;
 
Dim_sub2=cc; vector2=vec2; 
 
Entropy_PXP=zeros(Dim_sub,1);
 
for kk=1:Dim_sub
    
    kk
    rho0=zeros(Dim_sub2,Dim_sub2); 
  state=VV(:,kk);

for ii=1:Dim_sub2
    for jj=ii:Dim_sub2      

        
        find_index_ii=find(sum(abs(vector2(ii,:)-vector(:,1:num_site)),2)==0); %找出约化密度矩阵基矢对应的总基矢编号！
        find_index_jj=find(sum(abs(vector2(jj,:)-vector(:,1:num_site)),2)==0);
        
        [C,index_ii,index_jj]=intersect(vector(find_index_ii,num_site+1:end),vector(find_index_jj,num_site+1:end),'rows');
        
        
        
        rho0(ii,jj)=state(find_index_jj(index_jj),1)'*state(find_index_ii(index_ii),1); 
        rho0(jj,ii)=rho0(ii,jj)';
    end
end
Rho_temp=rho0;


Entropy_PXP(kk,1)=-trace(Rho_temp*logm(Rho_temp+0.000000000001));
warning('off')
end
 

Entropy_PXP=real(Entropy_PXP);
end
