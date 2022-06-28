function us=add_zeros(u,n)
L=length(u); z2a=floor((n-L)/2)+1;
if L<n
    if size(u,1)>1
        us=[zeros([1 z2a]) u' zeros([1 z2a-1])];
    else
        us=[zeros([1 z2a]) u zeros([1 z2a-1])];
    end
else
    if L>n
        if size(u,1)>1
            us=u(-z2a+1:L+z2a-1)';
        else
            us=u(-z2a+1:L+z2a-1);
        end
    else
        if L==n
            if size(u,1)>1
                us=u';
            else
                us=u;
            end
        end
    end
end
us=us(1:n);
    
