function [xshift,yshift] = circshift2centre_array(Xarray,Yarray)
% %% Do both rows seperately
x =sum(Xarray,2);
% x =max(Xarray,[],2);
p = find(x==max(x));
if rem(length(x),2)
    pshift = ((length(x)+1)/2)-p;
else
    pshift = ((length(x))/2)-p;
end
xshift = circshift(Xarray,pshift);

x =sum(Yarray,2);
% x =max(Yarray,[],2);
p = find(x==max(x));
if rem(length(x),2)
    pshift = ((length(x)+1)/2)-p;
else
    pshift = ((length(x))/2)-p;
end
yshift = circshift(Yarray,pshift);
%% Do both rows to off
% x =sum(Xarray,2);
% % % x =max(Xarray,[],2);
% 
% p = find(x==max(x));
% if rem(length(x),2)
%     pshift = ((length(x)+1)/2)-p;
% else
%     pshift = ((length(x))/2)-p;
% end
% xshift = circshift(Xarray,pshift);
% yshift = circshift(Yarray,pshift);
%% Do rows by the condition mean
% x =sum(Xarray+Yarray,2);
% % x =max(Xarray,[],2);
% p = find(x==max(x));
% if rem(length(x),2)
%     pshift = ((length(x)+1)/2)-p;
% else
%     pshift = ((length(x))/2)-p;
% end
% xshift = circshift(Xarray,pshift);
% yshift = circshift(Yarray,pshift);
%% Dont shift at all
% yshift = Yarray;
% xshift = Xarray
