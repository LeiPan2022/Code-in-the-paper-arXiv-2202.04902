function  vv_PXP=Find_v_PXP(L)

 [Dim_sub,vector,E_PXP,VV,find_index,index,Ham_PXP]=PXP_Ham(L);
 
 num_spinup=sum((vector==1),2);
 
 uu=1; 
for bb=1:1001
   vv=-(bb-1)/500;
   
   vv_b(bb,1)=vv;
    
  
BCS_wave=(uu.^(L-num_spinup).*vv.^(num_spinup));
BCS_wave=BCS_wave/norm(BCS_wave); 

Ener_Var1(bb)=BCS_wave'*Ham_PXP*BCS_wave;
 
Overlap1(bb,:)=BCS_wave(:,1)'*VV(:,1);
 
end

Ener_Var=real(Ener_Var1);

Ener_var_GS=min(Ener_Var1);
max(abs(Overlap1));

vv_PXP=vv_b(find(Ener_Var==min(Ener_Var1)));
 
end
 