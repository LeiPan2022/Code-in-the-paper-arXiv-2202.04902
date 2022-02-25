function  vv_PPXPP=Find_v_PPXPP(L)

 [Dim_sub_PPXPP,vector_PPXPP,E_PPXPP,VV_PPXPP,find_index_PPXPP,index_PPXPP,Ham_PPXPP]=PPXPP_Ham(L);
 
 num_spinup=sum((vector_PPXPP==1),2);
 
 uu=1; 
for bb=1:501
   vv=-(bb-1)/500;
   
   vv_b(bb,1)=vv;
    
  
BCS_wave=(uu.^(L-num_spinup).*vv.^(num_spinup));
BCS_wave=BCS_wave/norm(BCS_wave); 

 


Ener_Var1(bb)=BCS_wave'*Ham_PPXPP*BCS_wave;
 

Overlap1(bb,:)=BCS_wave(:,1)'*VV_PPXPP(:,1);
 
end

Ener_Var=real(Ener_Var1);

Ener_var_GS=min(Ener_Var1);
max(abs(Overlap1));

vv_PPXPP=vv_b(find(Ener_Var==min(Ener_Var1)))
 
end
 