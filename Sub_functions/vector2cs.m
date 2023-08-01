function cs=vector2cs(vector,N)

% Input data:
% vector: SHCs order in vetor as {S22,S21,C20,C21,C22;S33,S32,S31,C30,C31,C32,C33...}
% N: The maximum degree of SHCs

% Output data:
% cs: matrix |C\S| 

sz=0;
for iin=3:N+1
    num=2*(iin-1)+1;
    sc(iin,N+1-(num-1)/2:N+1+(num-1)/2)=vector(1+sz:num+sz);
    sz=sz+num;
end

[rows,cols] = size(sc);
lmax = rows -1;
c  = sc(:,lmax+1:2*lmax+1);
s  = [zeros(lmax+1,1) sc(:,1:lmax)];
cs = tril(c) + triu(rot90(s),1);

end