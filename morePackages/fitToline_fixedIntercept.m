function [b,Rsquared] = fitToline_fixedIntercept(xdata,ydata)
%fit to ydata = 1 + b*xdata
%   Detailed explanation goes here
fun = @(b,xdata) 1+b*xdata;

b = (sum(ydata.*xdata)- sum(xdata))/sum(xdata.^2);

%calculate the R square
ydata_bar = mean(ydata);
SS_tot = sum((ydata - repmat(ydata_bar,length(ydata),1)).^2);
SS_error = sum((ydata - fun(b,xdata)).^2);

Rsquared = 1 - SS_error / SS_tot;
end

