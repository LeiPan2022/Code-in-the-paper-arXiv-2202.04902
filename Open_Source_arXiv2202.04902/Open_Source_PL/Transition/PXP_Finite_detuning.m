function [vv_B,Gap_B,Overlap_B,Ener_B,Groundstate_B]=PXP_Finite_detuning(L,B)

  Dim=2^L;
 [Dim_sub,vector,E_PXP,VV,find_index,index,Ham_PXP,Sz_tot]=PXP_Ham_Sz(L);

 uu=1; 
for v_jj=1:161
    
    v_jj;
    vv=-0.4-(v_jj-1)/200;
    
  vv_b(v_jj)=vv;
 
phi_uv=cell(1,L+1);
 
phi_uv{1}=1;
phi_unit=[uu;vv]; phi_flip=[-vv;uu];

 

BCS_Basis=zeros(Dim_sub,L/2+1);BCS_Basis0=zeros(Dim,L/2+1);

for ll=1:L
    
    phi_uv{ll+1}=kron(phi_uv{ll},phi_unit);
    
end
    
BCS_Basis0(:,1)=phi_uv{L+1};

BCS_Basis(:,1)=BCS_Basis0(index,1);

BCS_Basis(:,1)=BCS_Basis(:,1)/norm(BCS_Basis(:,1)); 

nk=L/2; 
Momentum=exp(1i*2*pi*nk/L);
for ll=1:L
    sign=(Momentum)^(ll);
    
    BCS_Basis0(:,2)=BCS_Basis0(:,2)+sign*kron(kron(phi_uv{ll},phi_flip),phi_uv{L+1-ll});

end

BCS_Basis(:,2)=BCS_Basis0(index,2);
BCS_Basis(:,2)=BCS_Basis(:,2)/norm(BCS_Basis(:,2)); 

for ii=1:L-1
    for jj=ii+1:L
        sign=(Momentum)^(ii+jj);
        
        BCS_Basis0(:,3)=BCS_Basis0(:,3)+sign*kron(kron(kron(kron(phi_uv{ii},phi_flip),phi_uv{jj-ii}),phi_flip),phi_uv{L+1-jj});

    end
end

BCS_Basis(:,3)=BCS_Basis0(index,3);
BCS_Basis(:,3)=BCS_Basis(:,3)/norm(BCS_Basis(:,3)); 


if L>=6
for ii=1:L-2
    for jj=ii+1:L-1
        for kk=jj+1:L
            
        sign=(Momentum)^(ii+jj+kk);
        
        BCS_Basis0(:,4)=BCS_Basis0(:,4)+sign*kron(kron(kron(kron(kron(kron(phi_uv{ii},phi_flip),phi_uv{jj-ii}),phi_flip),phi_uv{kk-jj}),phi_flip),phi_uv{L+1-kk});
        
        end
    end
end

BCS_Basis(:,4)=BCS_Basis0(index,4);
BCS_Basis(:,4)=BCS_Basis(:,4)/norm(BCS_Basis(:,4)); 
end

if L>=8
for ii=1:L-3
    for jj=ii+1:L-2
        for kk=jj+1:L-1
            for ll=kk+1:L
        sign=(Momentum)^(ii+jj+kk+ll);
        
        BCS_Basis0(:,5)=BCS_Basis0(:,5)+sign*kron(kron(kron(kron(kron(kron(kron(kron(phi_uv{ii},phi_flip),phi_uv{jj-ii}),phi_flip),phi_uv{kk-jj}),phi_flip),phi_uv{ll-kk}),phi_flip),phi_uv{L+1-ll});
            
            end
        end
    end
end

BCS_Basis(:,5)=BCS_Basis0(index,5);
BCS_Basis(:,5)=BCS_Basis(:,5)/norm(BCS_Basis(:,5)); 
end


if L>=10
for ii=1:L-4
    for jj=ii+1:L-3
        for kk=jj+1:L-2
            for ll=kk+1:L-1
                for mm=ll+1:L
        sign=(Momentum)^(ii+jj+kk+ll+mm);
        
        BCS_Basis0(:,6)=BCS_Basis0(:,6)+sign*kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(phi_uv{ii},phi_flip),phi_uv{jj-ii}),phi_flip),phi_uv{kk-jj}),phi_flip),phi_uv{ll-kk}),phi_flip),phi_uv{mm-ll}),phi_flip),phi_uv{L+1-mm});
            
                end
            end
        end
    end
end

BCS_Basis(:,6)=BCS_Basis0(index,6);
BCS_Basis(:,6)=BCS_Basis(:,6)/norm(BCS_Basis(:,6)); 
end


if L>=12
for ii=1:L-5
    for jj=ii+1:L-4
        for kk=jj+1:L-3
            for ll=kk+1:L-2
                for mm=ll+1:L-1
                    for nn=mm+1:L
        sign=(Momentum)^(ii+jj+kk+ll+mm+nn);
        
        BCS_Basis0(:,7)=BCS_Basis0(:,7)+sign*kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(phi_uv{ii},phi_flip),phi_uv{jj-ii}),phi_flip),phi_uv{kk-jj}),phi_flip),phi_uv{ll-kk}),phi_flip),phi_uv{mm-ll}),phi_flip),phi_uv{nn-mm}),phi_flip),phi_uv{L+1-nn});
                    
                    end
                end
            end
        end
    end
end

BCS_Basis(:,7)=BCS_Basis0(index,7);
BCS_Basis(:,7)=BCS_Basis(:,7)/norm(BCS_Basis(:,7)); 
end


if L>=14
for ii=1:L-6
    for jj=ii+1:L-5
        for kk=jj+1:L-4
            for ll=kk+1:L-3
                for mm=ll+1:L-2
                    for nn=mm+1:L-1
                        for aa=nn+1:L
                            sign=(Momentum)^(ii+jj+kk+ll+mm+nn+aa);
                            BCS_Basis0(:,8)=BCS_Basis0(:,8)+sign*kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(phi_uv{ii},phi_flip),phi_uv{jj-ii}),phi_flip),phi_uv{kk-jj}),phi_flip),phi_uv{ll-kk}),phi_flip),phi_uv{mm-ll}),phi_flip),phi_uv{nn-mm}),phi_flip),phi_uv{aa-nn}),phi_flip),phi_uv{L+1-aa});
                        end
                    end
                end
            end
        end
    end
end

BCS_Basis(:,8)=BCS_Basis0(index,8);
BCS_Basis(:,8)=BCS_Basis(:,8)/norm(BCS_Basis(:,8));
end


if L>=16
for ii=1:L-7
    for jj=ii+1:L-6
        for kk=jj+1:L-5
            for ll=kk+1:L-4
                for mm=ll+1:L-3
                    for nn=mm+1:L-2
                        for aa=nn+1:L-1
                            for bb=aa+1:L
                                sign=(Momentum)^(ii+jj+kk+ll+mm+nn+aa+bb);
                                BCS_Basis0(:,9)=BCS_Basis0(:,9)+sign*kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(phi_uv{ii},phi_flip),phi_uv{jj-ii}),phi_flip),phi_uv{kk-jj}),phi_flip),phi_uv{ll-kk}),phi_flip),phi_uv{mm-ll}),phi_flip),phi_uv{nn-mm}),phi_flip),phi_uv{aa-nn}),phi_flip),phi_uv{bb-aa}),phi_flip),phi_uv{L+1-bb});
                            end
                        end
                    end
                end
            end
        end
    end
end

BCS_Basis(:,9)=BCS_Basis0(index,9);
BCS_Basis(:,9)=BCS_Basis(:,9)/norm(BCS_Basis(:,9));
end


if L>=18
for ii=1:L-8
    for jj=ii+1:L-7
        for kk=jj+1:L-6
            for ll=kk+1:L-5
                for mm=ll+1:L-4
                    for nn=mm+1:L-3
                        for aa=nn+1:L-2
                            for bb=aa+1:L-1
                                for cc=bb+1:L
                                    sign=(Momentum)^(ii+jj+kk+ll+mm+nn+aa+bb+cc);
                                    BCS_Basis0(:,10)=BCS_Basis0(:,10)+sign*kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(kron(phi_uv{ii},phi_flip),phi_uv{jj-ii}),phi_flip),phi_uv{kk-jj}),phi_flip),phi_uv{ll-kk}),phi_flip),phi_uv{mm-ll}),phi_flip),phi_uv{nn-mm}),phi_flip),phi_uv{aa-nn}),phi_flip),phi_uv{bb-aa}),phi_flip),phi_uv{cc-bb}),phi_flip),phi_uv{L+1-cc});
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
BCS_Basis(:,10)=BCS_Basis0(index,10);
BCS_Basis(:,10)=BCS_Basis(:,10)/norm(BCS_Basis(:,10));
end



BCS_Basis_New=zeros(Dim_sub,L/2+1);
BCS_Basis_New(:,1)=BCS_Basis(:,1);
for ii=2:L/2+1
for jj=1:ii-1
BCS_Basis_New(:,ii)=BCS_Basis_New(:,ii)-(BCS_Basis(:,ii)'*BCS_Basis_New(:,jj))/(BCS_Basis_New(:,jj)'*BCS_Basis_New(:,jj))*BCS_Basis_New(:,jj);
end
BCS_Basis_New(:,ii)=BCS_Basis_New(:,ii)+BCS_Basis(:,ii);
end

for k=1:L/2+1
BCS_Basis_New(:,k)=BCS_Basis_New(:,k)/norm(BCS_Basis_New(:,k));
end
 
 
 
 
 
 
 
 Ham_magnon2=zeros(L/2+1,L/2+1);
 
 for mm=1:L/2+1
 
    for nn=mm:L/2+1
    
        Ham_magnon2(mm,nn)=BCS_Basis_New(:,mm)'*(Ham_PXP-B*Sz_tot)*BCS_Basis_New(:,nn);
        Ham_magnon2(nn,mm)=Ham_magnon2(mm,nn)';

    end
end


[V_m2,Ener_m2]=eig(Ham_magnon2);
[DD_m2,J_m2]=sort(real(diag(Ener_m2)),'ascend');
% Ener_Magnons(:,cc)=diag(Ener_m2);
 Gap(v_jj,1)=DD_m2(2)-DD_m2(1);
 
 
 Ham_tot=Ham_PXP-B*Sz_tot;
[VV_B,Ener]=eig(Ham_tot);
%[VV,Ener]=eigs(Ham_tot,6,'sa');
[DDD,J]=sort(real(diag(Ener)),'ascend');
%Ener=sort(diag(Ener));
  
 
 New_vec=zeros(Dim_sub,L/2+1);

 for ii=1:L/2+1 
    for jj=1:L/2+1
        
    %New_vec(:,ii)=New_vec(:,ii)+V_m2(jj,J_m2(ii))*Vec_Eta_New(:,jj);
    New_vec(:,ii)=New_vec(:,ii)+V_m2(jj,J_m2(ii))*BCS_Basis_New(:,jj);
    
    end
    

 Overlap_ED_New(ii,v_jj)=New_vec(:,ii)'*VV_B(:,J(1));
  
 Ener_B(ii,v_jj)=real(New_vec(:,ii)'*(Ham_PXP-B*Sz_tot)*New_vec(:,ii));
 end
 
 
 
 
 GS_B(:,v_jj)=New_vec(:,1);      
  





end
 
 vv_B=vv_b(find(Ener_B(1,:)==min(Ener_B(1,:))))
 Overlap_B=Overlap_ED_New(1,(find(Ener_B(1,:)==min(Ener_B(1,:)))))
 Gap_B=Gap(find(Ener_B(1,:)==min(Ener_B(1,:))))
 Groundstate_B=GS_B(:,find(Ener_B(1,:)==min(Ener_B(1,:))));
 
end
 
 
 
 
 
 
 
 
 
 