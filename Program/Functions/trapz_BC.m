function output = trapz_BC(x,y)

if not(size(y,1)==size(x,1))
    y=y';
end

output=zeros(1,size(y,2));
for m=1:width(y)
    x2= x(~isnan(y(:,m)));
    y2= y(~isnan(y(:,m)),m);

    if x2(1)>x2(end)
        x2=flipud(x2);
        y2=flipud(y2);
    end

    output(:,m)=trapz(x2,linCorrect(y2));
end