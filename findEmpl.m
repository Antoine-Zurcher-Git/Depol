function [s] = findEmpl(num,xd,yd)
xt = table2array(num(:,1));
yt = table2array(num(:,2));

nF = -1;
for n=1:size(xt)
    if(xt(n) == xd && yt(n) == yd)
        nF = n;
    end
end
s = nF;

end