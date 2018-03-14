function xshift = circshift2centre(x,xrow)
        p = find(x(xrow,:)==max(x(xrow,:)));
        if rem(length(x),2)
            pshift = ((length(x)+1)/2)-p;
        else
            pshift = ((length(x))/2)-p;
        end
        xshift = circshift(x',pshift)';
% Do both rows seperately
% for i = 1:2
%     p = find(x(i,:)==max(x(i,:)));
%     pshift = (length(x)/2)-p;
%     xshift(i,:) = circshift(x(i,:)',pshift)';
% end