function [Dim_sub,vector,E_PXP,VV,find_index,index,Ham_PXP]=PXP_Ham(L)



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
vector=vec;
Ham_PXP=zeros(Dim_sub,Dim_sub); 
 


for ii=1:Dim_sub
    
     
    
    for kk=1:L
        
        
            vec1=zeros(1,L);
            vec1(1,kk)=(-1)^vector(ii,kk);
             
            find_state1=vector(ii,:)+vec1(1,:);
            
            jj=find(sum(abs(vector-find_state1),2)==0);
            
           if  ~isempty(jj)
            Ham_PXP(ii,jj)=1;
            Ham_PXP(jj,ii)=Ham_PXP(ii,jj)';
           end
    
    
    end
               
               
end
                           
       

[VV,Ener]=eig(Ham_PXP);
[DDD,J]=sort(diag(Ener),'ascend');
E_PXP=diag(Ener);
%VV=VV(:,J);

end



