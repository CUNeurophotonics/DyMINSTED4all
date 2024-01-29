function F = D2LorentianFunction(x,xdata)


%  F = x(1)*pi*(1/2*x(3)./((xdata(1,:,1).^2 -x(2)).^2 + (1/2*x(3))^2))+x(6);
 
%  F = x(1)*pi/2*(x(3)+x(5))/2*(1/2*x(3)./((xdata(:,:,1)-x(2)).^2 + (1/2*x(3))^2))...
%      .*(1/2*x(5)./((xdata(:,:,2)-x(4)).^2 + (1/2*x(5))^2))+x(6);
 
 F = x(1)*(x(3)*x(5))/4*(1/2*x(3)./((xdata(:,:,1)-x(2)).^2 + (1/2*x(3))^2))...
     .*(1/2*x(5)./((xdata(:,:,2)-x(4)).^2 + (1/2*x(5))^2))+x(6);
 
 
 % + 1/2*x(5)./((xdata(:,:,2)-x(4)).^2 + (1/2*x(5))^2) 
 
 %x(1) is amplidute
 %x(2) is offset in the x direction
 %x(3) is FWHM along x direction
 %x(4) is offset in the y direction
 %x(5) is FWHM along y direction
 %x(6) is offset
 
 %xdata(:,:,1) are values of x coordinates
 %xdata(:,:,2) are values of y coordinates