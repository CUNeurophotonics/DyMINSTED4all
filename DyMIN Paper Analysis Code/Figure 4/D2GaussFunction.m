function F = D2GaussFunction(x,xdata)
 F = x(1)*exp(   -((xdata(:,:,1)-x(2)).^2/(2*x(3)^2) + (xdata(:,:,2)-x(4)).^2/(2*x(5)^2) )    )+x(6);
 
 %x(1) is amplidute
 %x(2) is offset in the x direction
 %x(3) is waist along x direction
 %x(4) is offset in the y direction
 %x(5) is waist along y direction
 %x(6) is offset
 
 %xdata(:,:,1) are values of x coordinates
 %xdata(:,:,2) are values of y coordinates