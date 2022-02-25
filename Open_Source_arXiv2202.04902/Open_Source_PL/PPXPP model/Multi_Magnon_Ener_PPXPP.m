function [Max_Overlap2,Max_index2,Ener_Multi_Magnon_PPXPP,Sz_average_Magnons]=Multi_Magnon_Ener_PPXPP(L)

uu=1;  
Dim=2^L;
[Dim_sub_PPXPP,vector_PPXPP,E_PPXPP,VV_PPXPP,find_index_PPXPP,index_PPXPP,Ham_PPXPP]=PPXPP_Ham(L);
[Ener_k,VV_k,Sz_average_Mom]=PPXPP_Ener_Sz_Momentum(L);


vv=Find_v_PPXPP(L);
phi_uv=cell(1,L+1);
 
phi_uv{1}=1;
phi_unit=[uu;vv]; phi_flip=[-vv;uu];

BCS_Basis=zeros(Dim_sub_PPXPP,1+L); BCS_Basis0=zeros(Dim,1+L);



for ll=1:L
    
    phi_uv{ll+1}=kron(phi_uv{ll},phi_unit);
    
end
    
BCS_Basis0(:,1)=phi_uv{L+1};

BCS_Basis(:,1)=BCS_Basis0(index_PPXPP,1);

BCS_Basis(:,1)=BCS_Basis(:,1)/norm(BCS_Basis(:,1)); 


Momentum1=exp(1i*2*pi/3);

Momentum2=exp(1i*pi*4/3);

for ll=1:L
    sign=(Momentum1)^(ll);
    
    BCS_Basis0(:,2)=BCS_Basis0(:,2)+sign*kron(kron(phi_uv{ll},phi_flip),phi_uv{L+1-ll});

end

BCS_Basis(:,2)=BCS_Basis0(index_PPXPP,2);
BCS_Basis(:,2)=BCS_Basis(:,2)/norm(BCS_Basis(:,2)); 


for ll=1:L
    sign=(Momentum2)^(ll);
    
    BCS_Basis0(:,3)=BCS_Basis0(:,3)+sign*kron(kron(phi_uv{ll},phi_flip),phi_uv{L+1-ll});

end

BCS_Basis(:,3)=BCS_Basis0(index_PPXPP,3);
BCS_Basis(:,3)=BCS_Basis(:,3)/norm(BCS_Basis(:,3)); 


for ii=1:L-1
    for jj=ii+1:L
                  
        sign=(Momentum1)^(ii+jj);
        
        BCS_Basis0(:,4)=BCS_Basis0(:,4)+sign*kron(kron(kron(kron(phi_uv{ii},phi_flip),phi_uv{jj-ii}),phi_flip),phi_uv{L+1-jj});
                 
    end
end

BCS_Basis(:,4)=BCS_Basis0(index_PPXPP,4);
BCS_Basis(:,4)=BCS_Basis(:,4)/norm(BCS_Basis(:,4)); 


for ii=1:L-1
    for jj=ii+1:L
                  
        sign=(Momentum2)^(ii+jj);
        
        BCS_Basis0(:,5)=BCS_Basis0(:,5)+sign*kron(kron(kron(kron(phi_uv{ii},phi_flip),phi_uv{jj-ii}),phi_flip),phi_uv{L+1-jj});
                 
    end
end

BCS_Basis(:,5)=BCS_Basis0(index_PPXPP,5);
BCS_Basis(:,5)=BCS_Basis(:,5)/norm(BCS_Basis(:,5)); 

for ii=1:L-1
    for jj=ii+1:L
                  
        sign=(Momentum1)^ii*(Momentum2)^jj+(Momentum1)^jj*(Momentum2)^ii;
        
        BCS_Basis0(:,6)=BCS_Basis0(:,6)+sign*kron(kron(kron(kron(phi_uv{ii},phi_flip),phi_uv{jj-ii}),phi_flip),phi_uv{L+1-jj});
                 
    end
end

BCS_Basis(:,6)=BCS_Basis0(index_PPXPP,6);
BCS_Basis(:,6)=BCS_Basis(:,6)/norm(BCS_Basis(:,6)); 

 
for ii=1:L-2
    for jj=ii+1:L-1
        for kk=jj+1:L
            
        sign=(Momentum2)^ii*(Momentum1)^(jj+kk)+(Momentum2)^jj*(Momentum1)^(ii+kk)+(Momentum2)^kk*(Momentum1)^(ii+jj);
        
        BCS_Basis0(:,7)=BCS_Basis0(:,7)+sign*kron(kron(kron(kron(kron(kron(phi_uv{ii},phi_flip),phi_uv{jj-ii}),phi_flip),phi_uv{kk-jj}),phi_flip),phi_uv{L+1-kk});
        
        end
    end
end

BCS_Basis(:,7)=BCS_Basis0(index_PPXPP,7);
BCS_Basis(:,7)=BCS_Basis(:,7)/norm(BCS_Basis(:,7)); 

 
for ii=1:L-2
    for jj=ii+1:L-1
        for kk=jj+1:L
            
        sign=(Momentum1)^ii*(Momentum2)^(jj+kk)+(Momentum1)^jj*(Momentum2)^(ii+kk)+(Momentum1)^kk*(Momentum2)^(ii+jj);
        
        BCS_Basis0(:,8)=BCS_Basis0(:,8)+sign*kron(kron(kron(kron(kron(kron(phi_uv{ii},phi_flip),phi_uv{jj-ii}),phi_flip),phi_uv{kk-jj}),phi_flip),phi_uv{L+1-kk});
        
        end
    end
end

BCS_Basis(:,8)=BCS_Basis0(index_PPXPP,8);
BCS_Basis(:,8)=BCS_Basis(:,8)/norm(BCS_Basis(:,8)); 


for ii=1:L-2
    for jj=ii+1:L-1
        for kk=jj+1:L
            
        sign=(Momentum1)^(ii+jj+kk);
        
        BCS_Basis0(:,9)=BCS_Basis0(:,9)+sign*kron(kron(kron(kron(kron(kron(phi_uv{ii},phi_flip),phi_uv{jj-ii}),phi_flip),phi_uv{kk-jj}),phi_flip),phi_uv{L+1-kk});
        
        end
    end
end

BCS_Basis(:,9)=BCS_Basis0(index_PPXPP,9);
BCS_Basis(:,9)=BCS_Basis(:,9)/norm(BCS_Basis(:,9)); 



for ii=1:L-2
    for jj=ii+1:L-1
        for kk=jj+1:L
            
        sign=(Momentum2)^(ii+jj+kk);
        
        BCS_Basis0(:,10)=BCS_Basis0(:,10)+sign*kron(kron(kron(kron(kron(kron(phi_uv{ii},phi_flip),phi_uv{jj-ii}),phi_flip),phi_uv{kk-jj}),phi_flip),phi_uv{L+1-kk});
        
        end
    end
end

BCS_Basis(:,10)=BCS_Basis0(index_PPXPP,10);
BCS_Basis(:,10)=BCS_Basis(:,10)/norm(BCS_Basis(:,10)); 


for ii=1:L-3
    for jj=ii+1:L-2
        for kk=jj+1:L-1
            for ll=kk+1:L
  sign=(Momentum1)^(ii+jj)*(Momentum2)^(kk+ll)+(Momentum1)^(ii+kk)*(Momentum2)^(jj+ll)+(Momentum1)^(ii+ll)*(Momentum2)^(jj+kk)...
      +(Momentum1)^(jj+kk)*(Momentum2)^(ii+ll)+(Momentum1)^(jj+ll)*(Momentum2)^(ii+kk)+(Momentum1)^(kk+ll)*(Momentum2)^(ii+jj);

        
        BCS_Basis0(:,11)=BCS_Basis0(:,11)+sign*kron(kron(kron(kron(kron(kron(kron(kron(phi_uv{ii},phi_flip),phi_uv{jj-ii}),phi_flip),phi_uv{kk-jj}),phi_flip),phi_uv{ll-kk}),phi_flip),phi_uv{L+1-ll});
            
            end
        end
    end
end

BCS_Basis(:,11)=BCS_Basis0(index_PPXPP,11);
BCS_Basis(:,11)=BCS_Basis(:,11)/norm(BCS_Basis(:,11)); 



for ii=1:L-3
    for jj=ii+1:L-2
        for kk=jj+1:L-1
            for ll=kk+1:L
  sign=(Momentum2)^ii*(Momentum1)^(jj+kk+ll)+(Momentum2)^jj*(Momentum1)^(ii+kk+ll)+(Momentum2)^kk*(Momentum1)^(ii+jj+ll)+(Momentum2)^ll*(Momentum1)^(ii+jj+kk);

        
        BCS_Basis0(:,12)=BCS_Basis0(:,12)+sign*kron(kron(kron(kron(kron(kron(kron(kron(phi_uv{ii},phi_flip),phi_uv{jj-ii}),phi_flip),phi_uv{kk-jj}),phi_flip),phi_uv{ll-kk}),phi_flip),phi_uv{L+1-ll});
            
            end
        end
    end
end

BCS_Basis(:,12)=BCS_Basis0(index_PPXPP,12);
BCS_Basis(:,12)=BCS_Basis(:,12)/norm(BCS_Basis(:,12)); 


for ii=1:L-3
    for jj=ii+1:L-2
        for kk=jj+1:L-1
            for ll=kk+1:L
  sign=(Momentum1)^ii*(Momentum2)^(jj+kk+ll)+(Momentum1)^jj*(Momentum2)^(ii+kk+ll)+(Momentum1)^kk*(Momentum2)^(ii+jj+ll)+(Momentum1)^ll*(Momentum2)^(ii+jj+kk);

        
        BCS_Basis0(:,13)=BCS_Basis0(:,13)+sign*kron(kron(kron(kron(kron(kron(kron(kron(phi_uv{ii},phi_flip),phi_uv{jj-ii}),phi_flip),phi_uv{kk-jj}),phi_flip),phi_uv{ll-kk}),phi_flip),phi_uv{L+1-ll});
            
            end
        end
    end
end

BCS_Basis(:,13)=BCS_Basis0(index_PPXPP,13);
BCS_Basis(:,13)=BCS_Basis(:,13)/norm(BCS_Basis(:,13)); 


for ii=1:L-3
    for jj=ii+1:L-2
        for kk=jj+1:L-1
            for ll=kk+1:L
  sign=(Momentum1)^(ii+jj+kk+ll);
        
        BCS_Basis0(:,14)=BCS_Basis0(:,14)+sign*kron(kron(kron(kron(kron(kron(kron(kron(phi_uv{ii},phi_flip),phi_uv{jj-ii}),phi_flip),phi_uv{kk-jj}),phi_flip),phi_uv{ll-kk}),phi_flip),phi_uv{L+1-ll});
            
            end
        end
    end
end

BCS_Basis(:,14)=BCS_Basis0(index_PPXPP,14);
BCS_Basis(:,14)=BCS_Basis(:,14)/norm(BCS_Basis(:,14)); 


for ii=1:L-3
    for jj=ii+1:L-2
        for kk=jj+1:L-1
            for ll=kk+1:L
  sign=(Momentum2)^(ii+jj+kk+ll);
        
        BCS_Basis0(:,15)=BCS_Basis0(:,15)+sign*kron(kron(kron(kron(kron(kron(kron(kron(phi_uv{ii},phi_flip),phi_uv{jj-ii}),phi_flip),phi_uv{kk-jj}),phi_flip),phi_uv{ll-kk}),phi_flip),phi_uv{L+1-ll});
            
            end
        end
    end
end

BCS_Basis(:,15)=BCS_Basis0(index_PPXPP,15);
BCS_Basis(:,15)=BCS_Basis(:,15)/norm(BCS_Basis(:,15)); 



BCS_Basis_New=zeros(Dim_sub_PPXPP,15);
BCS_Basis_New(:,1)=BCS_Basis(:,1);
for ii=2:15
for jj=1:ii-1
BCS_Basis_New(:,ii)=BCS_Basis_New(:,ii)-(BCS_Basis(:,ii)'*BCS_Basis_New(:,jj))/(BCS_Basis_New(:,jj)'*BCS_Basis_New(:,jj))*BCS_Basis_New(:,jj);
end
BCS_Basis_New(:,ii)=BCS_Basis_New(:,ii)+BCS_Basis(:,ii);
end

for k=1:15
BCS_Basis_New(:,k)=BCS_Basis_New(:,k)/norm(BCS_Basis_New(:,k));
end
 
  
for mm=1:15
    mm
    for nn=1:15
    
        Ham_Magnon2(mm,nn)=BCS_Basis_New(:,mm)'*Ham_PPXPP*BCS_Basis_New(:,nn);
         
    end
end


[V_M2,Ener_M2]=eig(Ham_Magnon2);
[DD_M2,J_M2]=sort(real(diag(Ener_M2)),'ascend');
 DD_M2;
V_M2(:,J_M2(1))



New_vec=zeros(Dim_sub_PPXPP,15);
for ii=1:15
    
    for jj=1:15
        
    New_vec(:,ii)=New_vec(:,ii)+V_M2(jj,J_M2(ii))*BCS_Basis_New(:,jj);
    
    end
    
end


for ii=1:15

     New_Overlap(ii,:)=abs(New_vec(:,ii)'*VV_k);
     Max_Overlap(ii,1)=max(New_Overlap(ii,:));
     
     Max_index(ii,1)=find(New_Overlap(ii,:)==max(New_Overlap(ii,:)));
     
     New_Overlap2(ii,:)=abs(BCS_Basis_New(:,ii)'*VV_k);
    
     Max_Overlap2(ii,1)=max(New_Overlap2(ii,:));
     
     Max_index2(ii,1)=find(New_Overlap2(ii,:)==max(New_Overlap2(ii,:)));
end

Max_Overlap;
Max_index;
Max_Overlap2;
Max_index2;

Sz_1average_Magnons=zeros(L,15);
for jj=1:15
 
    Ener_Multi_Magnon_PPXPP(1,jj)=BCS_Basis_New(:,jj)'*Ham_PPXPP*BCS_Basis_New(:,jj);
    for kk=1:L
   
        Sz_1average_Magnons(kk,jj)=(2*vector_PPXPP(:,kk)-1)'*(BCS_Basis_New(:,jj).*conj(BCS_Basis_New(:,jj)));
        
    end
end

Sz_average_Magnons(1,:)=sum(Sz_1average_Magnons)/L;

end

 
