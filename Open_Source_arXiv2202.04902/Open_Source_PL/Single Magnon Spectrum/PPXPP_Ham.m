function [Dim_sub_PPXPP,vector_PPXPP,E_PPXPP,VV_PPXPP,find_index_PPXPP,index_PPXPP,Ham_PPXPP]=PPXPP_Ham(L)



Dim=2^L;

temp=zeros(1,L); temp1=zeros(1,L); temp2=zeros(1,L);  
 bb=0;
for jj=1:Dim
    
obj=dec2base(jj-1,2,L);  
temp=(obj=='1')*1+(obj=='0')*0;

temp1(2:end)=temp(1:end-1);
temp1(1)=temp(end);
temp2(3:end)=temp(1:end-2); 
temp2(2)=temp(end); temp2(1)=temp(end-1);

 

kk1=sum(temp.*temp1);  kk2=sum(temp.*temp2);  
   if kk1==0 && kk2==0  
       bb=bb+1;
      vec_PPXPP(bb,:)=temp(1,:); 
      Sz_temp(bb)=sum(vec_PPXPP(bb,:)-1/2);
      index_PPXPP(bb,:)=jj;
      find_index_PPXPP(jj,1)=bb;
   end
end

 

Dim_sub_PPXPP=bb;
vector_PPXPP=vec_PPXPP;
Ham_PPXPP=zeros(Dim_sub_PPXPP,Dim_sub_PPXPP); 
 


for ii=1:Dim_sub_PPXPP
    
     
    
    for kk=1:L
        
        
            vec1=zeros(1,L);
            vec1(1,kk)=(-1)^vector_PPXPP(ii,kk);          
            find_state1=vector_PPXPP(ii,:)+vec1(1,:);
            
            jj=find(sum(abs(vector_PPXPP-find_state1),2)==0);
            
           if  ~isempty(jj)
            Ham_PPXPP(ii,jj)=1;
            Ham_PPXPP(jj,ii)=Ham_PPXPP(ii,jj)';
           end
    
    
    end
               
               
end                           
      

[VV_PPXPP,Ener_PPXPP]=eig(Ham_PPXPP);
[DDD,J]=sort(diag(Ener_PPXPP),'ascend');
E_PPXPP=diag(Ener_PPXPP);
%VV=VV(:,J);

end



