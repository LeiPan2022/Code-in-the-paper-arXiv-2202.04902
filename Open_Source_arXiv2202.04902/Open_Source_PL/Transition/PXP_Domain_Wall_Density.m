function [Domain_wall_mean,Domain_wall_flu]=PXP_Domain_Wall_Density(L,Groundstate_B)

[Dim_sub,vector,E_PXP,VV,find_index,index,Ham_PXP,Sz_tot]=PXP_Ham_Sz(L);

for kk=1:Dim_sub
    
    Domain_num(kk,1)=sum([diff(vector(kk,:)),vector(kk,1)-vector(kk,end)]==0);

end
   
Domain_wall_mean=abs(Groundstate_B').^2*Domain_num/L;

Domain_wall_flu=(abs(Groundstate_B').^2*(Domain_num/L).^2-Domain_wall_mean^2)/1;
   
end

