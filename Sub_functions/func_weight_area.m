function w=func_weight_area(th1,th2,d)

% Input data:
% th1, th2: Co-latitude 

% Output data:
% w: The weight of point on the sphere

n=180/d;
L_n=(0:n/2-1); % n should be an even number

w=zeros(round((th2-th1)/d+1),1);
j=0;
for th=th1:d:th2
    j=j+1;
    w(j)=4*pi/n^2*sin(th*pi/180)*sum(sin((2*L_n+1)*th*pi/180)./(2*L_n+1));
end

end

