%%calcline
%%arr 11/5/2007

function a1 = calcline(b,a1,b1,a2,b2)

if b2==b1;
    a1=a1;
else
    a1=a1+round(double(b-b1)*double(a2-a1)/double(b2-b1));
end
