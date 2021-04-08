%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FINITE DIFFERENCE DISCRETIZATION OF THE TWO DIMENSIONAL LAPLACE EQUATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISCRETIZATION: : SECOND ORDER CENTRAL DIFFERENCE SCHEME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE WRITTEN BY: NAHOM ALEMSEGED WORKU
% GENERAL EQUATION: D2P/DX2 + D2P/DY2 = 0
% BOUNDARY CONDITION: P = 0 @ X=0, P = Y+1 @ X = 5, DP/DY = 0 @Y = 0,Y = 5
% DOMAIN = [0 X 5] [0 X 5]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLE DECLERATION
Xp = 5;    Yp = 5;
nx = 6;    ny = 6;
dx = Xp/(nx-1);    dy = Yp/(ny - 1);
P = zeros(nx,ny);
Pd = zeros(nx,ny);
X = [0:dx:Xp];          Y = [0:dy:Yp]';
maxit = 500;
% BOUNDARY CONDITIONS
YPP = [Yp-1:-dy:1];
P([2:5],end)  = YPP + 1;
    for j = Yp+1:-dy:2
       P(1,j-1) = P(1,j);
       P(end,j-1) = P(end,j)
    end
%    FIND PRESSURE TERMS
for it = 1:maxit
    pd = P;
    for i=2:nx-1
        for j=2:ny-1
            P(i,j) = (((dy^2).*(Pd(i-1,j) + Pd(i+1,j))) +  ((dx^2).*(Pd(i,j-1) + Pd(i,j+1)))) / (2.*((dx^2) + (dy^2)))
        end
    end
end

 contourf(X,Y,P)
 colorbar;
grid on;
% quiver(X,Y,P)




    
    
    
    
