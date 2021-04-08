%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FINITE DIFFERENCE DISCRETIZATION OF THE TWO DIMENSIONAL POISSON'S  EQUATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISCRETIZATION: : SECOND ORDER CENTRAL DIFFERENCE SCHEME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE WRITTEN BY: NAHOM ALEMSEGED WORKU
% GENERAL EQUATION: D2P/DX2 + D2P/DY2 = b
% BOUNDARY CONDITION: P = 0 @ X=0, P = 0 @ X = 5, DP/DY = 0 @Y = 0,Y = 5
% SOURCE TERMS: b = 100 @ X = X/4 AND Y = Y/4
%                          b = -100 @ X = 3/4 AND Y = Y/4
%                          b = 0, EVERYWHERE ELSE
% DOMAIN = [0 X 5] [0 X 5]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLE DECLERATION
Xp = 8;    Yp = 8;
nx = 9;    ny = 9;
dx = Xp/(nx-1);    dy = Yp/(ny - 1);
P = zeros(nx,ny);
b = zeros(nx,ny);
X = [0:dx:Xp];          Y = [0:dy:Yp]';
maxit = 500;
% BOUNDARY CONDITIONS
b(Xp/4,Yp/4) = 100;
b(3*Xp/4,3*Yp/4) = -100;

%    FIND PRESSURE TERMS
for it = 1:maxit
    for i=2:nx-1
        for j=2:ny-1
                 P(i,j) = (((dy^2).*(P(i-1,j) + P(i+1,j))) +  ((dx^2).*(P(i,j-1) + P(i,j+1))) - (b(i,j) * (dx^2) * (dy^2))) / (2.*((dx^2) + (dy^2)))
        end
    end
end

 contourf(X,Y,P)
 colorbar;
grid on;
% quiver(X,Y,P)




    
    
    
    
