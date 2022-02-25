function [Ener_k,VV_k,Sz_average_Mom]=PPXPP_Ener_Sz_Momentum(L)


[Dim_sub_PPXPP,vector_PPXPP,E_PPXPP,VV_PPXPP,find_index_PPXPP,index_PPXPP,Ham_PPXPP]=PPXPP_Ham(L);

T_Matrix=zeros(Dim_sub_PPXPP,Dim_sub_PPXPP);
T_Matrix2=zeros(Dim_sub_PPXPP,Dim_sub_PPXPP);

 T_Matrix_set=cell(L,1);
    
        
        
        for ii=1:L
            ii;
            T_Matrix_temp=sparse(Dim_sub_PPXPP,Dim_sub_PPXPP);
            
            for jj=1:Dim_sub_PPXPP
                
        jj2=find_index_PPXPP(basis_to_index(circshift(vector_PPXPP(jj,:),ii-1))+1);  % ii=1，表示平移零次
        
             T_Matrix_temp(jj2,jj)=1;
             
            end
        
         T_Matrix_set{ii,1}=T_Matrix_temp;

        end

% Momentum subspace

U_matrix_K=cell(L,1); % Momentum transformation matrix

for kk=1:L  
    tic
    kk;
    U_matrix_K{kk,1}=sparse(Dim_sub_PPXPP,Dim_sub_PPXPP);
  for jj=1:L 
    
    U_matrix_K{kk,1}=U_matrix_K{kk,1}+exp(1i*2*pi*(kk-1)/L*(jj-1))*T_Matrix_set{jj,1};
    
  end
  
  if kk==1
      U_matrix=orth(full(U_matrix_K{kk,1}));
      Dim_block(kk,1)=rank(U_matrix);
  end
  
  if kk>1
      U_matrix=[U_matrix,orth(full(U_matrix_K{kk,1}))];
      Dim_block(kk,1)=rank(orth(full(U_matrix_K{kk,1})));
  end
  toc
end
    Dim_block;
    
    Ham_PPXPP2=U_matrix'*Ham_PPXPP*U_matrix;
    Ham_PPXPP2=(Ham_PPXPP2+Ham_PPXPP2')/2;
    
    
%     VV_k=cell(L+1/2,1);
%     Ener_k=cell(L+1/2,1);
    

    
%for kk=1:L/2+1
for kk=1:L
       kk;
    if kk==1
        [VV2,Ener2]=eig(Ham_PPXPP2(1:Dim_block(1),1:Dim_block(1)));
        [DDD2,J2]=sort(diag(Ener2),'ascend');
    
        VV_k(:,1:Dim_block(1))=U_matrix(:,1:Dim_block(1))*VV2;
        Ener_k(1:Dim_block(1),1)=diag(Ener2);
    end
    
    
    if kk>1
        [VV3,Ener3]=eig(Ham_PPXPP2(sum(Dim_block(1:kk-1))+1:sum(Dim_block(1:kk-1))+Dim_block(kk),sum(Dim_block(1:kk-1))+1:sum(Dim_block(1:kk-1))+Dim_block(kk)));
        [DDD3,J3]=sort(diag(Ener3),'ascend');
        
        VV_k(:,sum(Dim_block(1:kk-1))+1:sum(Dim_block(1:kk-1))+Dim_block(kk))=U_matrix(:,sum(Dim_block(1:kk-1))+1:sum(Dim_block(1:kk-1))+Dim_block(kk))*VV3;
        Ener_k(sum(Dim_block(1:kk-1))+1:sum(Dim_block(1:kk-1))+Dim_block(kk),1)=diag(Ener3);
    end
    
    
 
    
end   




Sz_1average_Mom=zeros(L,Dim_sub_PPXPP);
for jj=1:Dim_sub_PPXPP
 
    for kk=1:L
   
        Sz_1average_Mom(kk,jj)=(2*vector_PPXPP(:,kk)-1)'*(VV_k(:,jj).*conj(VV_k(:,jj)));
        
    end
end

Sz_average_Mom(1,:)=sum(Sz_1average_Mom)/L;
