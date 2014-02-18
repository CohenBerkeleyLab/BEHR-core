%stdSCIA.m
%Takes the seasonal data and gets a n_std 

function k1 = stdSCIA(a1, b1, c1, d1, e1, f1)
z=[];
n=240;            %number of rows
m=440;            %number of columns

for j=1:n;          
    t1=[];
    for i=1:m,      
        R=[a1(j,i), b1(j,i), c1(j,i), d1(j,i), e1(j,i), f1(j,i)];
        t=n_std(R);
        t1=[t1 t];
    end
    z=[z;t1];
end

k1=z;  