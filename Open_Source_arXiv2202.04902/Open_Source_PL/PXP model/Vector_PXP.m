function [Dim_sub,vector_PXP,Dim_sub2,vector_PXP2]=Vector_PXP(L)


Dim=2^L;

temp=zeros(1,L); temp1=zeros(1,L); temp2=zeros(1,L); temp3=zeros(1,L);
 bb=0;
for jj=1:Dim
    jj;
obj=dec2base(jj-1,2,L);  
temp=(obj=='1')*1+(obj=='0')*0;

temp1(2:end)=temp(1:end-1);
temp1(1)=temp(end);

 

kk1=sum(temp.*temp1);  
   if kk1==0  
       bb=bb+1;
      vec(bb,:)=temp(1,:); 
      Sz_temp(bb)=sum(vec(bb,:)-1/2);
      index(bb,:)=jj;
      find_index(jj,1)=bb;
   end
end

 

Dim_sub=bb;
vector_PXP=vec;


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
 
Dim_sub2=cc; vector_PXP2=vec2; 


