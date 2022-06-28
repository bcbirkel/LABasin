function [ys,indx]=sort_mex(x,y,var,res,n,m)

cnt=zeros(m,1);
gamma_sum=zeros(m,1);
var_sum=zeros(m,1);
var2_sum=zeros(m,1);
gammah=zeros(m,1);
c0=zeros(m,1);
if n < 1000
  n_cnt=max(1,round(0.05*n));
else
  n_cnt=max(1,round(0.01*n));
end
head_indx=zeros(n,m);

i0=1;
for i=1:n
	xi=x(i+1:n);yi=y(i+1:n);var_i=var(i+1:n);ni=n-i;
   d=sqrt( (x(i)-xi).^2+(y(i)-yi).^2);
   [dsort,indx_sort]=sort(d);
   dsort=dsort(i0:ni);								% excluding
   indx_sort=indx_sort(i0:ni);
   var_sort=var_i(indx_sort);       
   dindx=floor(dsort/res)+1;
   k_acc=0;								% accumulated k index
   for k=1:m
     indxk=find(dindx(k_acc+1:ni+1-i0) == k) + k_acc ;
     nk=length(indxk);
     if nk > 0
		% tail portion
        cnt(k)=cnt(k)+nk;
        k_acc=k_acc+nk;
        var_dif= var(i) - var_sort(indxk);	
        var_sum(k)=var_sum(k)+sum(var_sort(indxk));
        var2_sum(k)=var2_sum(k)+var_sort(indxk)'*var_sort(indxk);
        gamma_sum(k)=gamma_sum(k)+ var_dif'*var_dif;
		% head point index
		  head_indx(i,k)=1;
      end
   end
end

ys=y;
indx=head_indx;