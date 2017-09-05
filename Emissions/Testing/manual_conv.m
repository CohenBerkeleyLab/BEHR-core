function [  ] = manual_conv( u, v )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


m = length(u);
n = length(v);
w = nan(1, m+n-1);

for k=1:length(w)
    j = max(1, k+1-n):min(k,m);
    j2 = k - j + 1;
    f=figure; plot(u(j)/max(u(j))); hold on; plot(v(j2)/max(v(j2)));
    title(sprintf('k = %d',k))
    input('enter');
    close(f);
end

end

