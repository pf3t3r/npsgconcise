function xout = modelEquation(F,t0,h,tFinal,x0)

x = x0;
xout = x;

for t= t0 : h : tFinal-h
    s1 = F(t,x);
    s2 = F(t+h/2,x+h*s1/2);
    s3 = F(t+h/2,x+h*s2/2);
    s4 = F(t+h,x+h*s3);
    x = x + h*(s1 + 2*s2 + 2*s3 + s4)/6;
    xout = [xout; x];
end
