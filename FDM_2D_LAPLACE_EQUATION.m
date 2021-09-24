%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FINITE DIFFERENCE DISCRETIZATION OF THE TWO DIMENSIONAL LAPLACE EQUATION
%EXPLICIT METHOD DEMO PROGRAM.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISCRETIZATION: : SECOND ORDER CENTRAL DIFFERENCE SCHEME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE WRITTEN BY: NAHOM ALEMSEGED WORKU
% GENERAL EQUATION: D2P/DX2 + D2P/DY2 = 0
% BOUNDARY CONDITION: P' = 0 @ X=0, P' = 0 @ X = L, P' = 0 @Y = 0,Y = L
% SIMPLE LAPLACE EQUATION WITH A BARRIER AT THE CENTER WHICH CAN ALXO BE ADOPTED FOR SEEPAGE ANALYSIS OF DAMS
% UPSTREAM PRESSURE IS 44.9 AND DOWNSTREAM HYDRAULIC HEAD IS 36.5
% DOMAIN = [0 X Xp] [0 X Yp]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLE DECLERATION
Xp = 1;    Yp = 1;
nx = 50;    ny = 50;
dx = Xp/(nx-1);    dy = Yp/(ny - 1);
P = zeros(nx,ny);
pd = zeros(nx,ny);
X = [0:dx:Xp];          Y = [0:dy:Yp]';
maxit = 10000; % MAXIMUM ITERATION
% BOUNDARY CONDITIONS

P(end,:)  = 44.9;

    %    FIND PRESSURE TERMS
for it = 1:maxit
%     pd = P;
    for i=2:nx-1
        for j=2:ny-1
            pd = P;
%             P(i,j) = (((dy.^2).*(pd(i-1,j) + pd(i+1,j))) +  ((dx.^2).*(pd(i,j-1) + pd(i,j+1)))) / (2.*((dx.^2) + (dy.^2)));
              P(i,j) = (((dy.^2).*(P(i-1,j) + P(i+1,j))) +  ((dx.^2).*(P(i,j-1) + P(i,j+1)))) / (2.*((dx.^2) + (dy.^2)));
            P(1,:) = P(2,:);
            P(:,1) = P(:,2);
            P(:,end) = P(:,end-1);
            P(end,[1:25])  = 44.9;
            P(end,[26:end])  = 36.5;
       
            P(end,[20:25]) =  P(end-1,[20:25]);
%             P([end-5, end],20) = 0;
%             P([end-5, end],21) = P([end-5, end],19);
            P([end/2:end],[19:21]) = 0;
%             P([end/2, end],21) = P([end-5, end],19);

        end
    end
end
%   if i<25 && j==50
%                P(end,:)  = 44.9;
%              elseif i >= 25 && j == 50
%                   P(end,:)  = 35.6;
%            end
[c,h] = contourf(X,Y,P);
set(h, 'edgecolor','none');
colormap('jet')
colorbar;
grid on;

% quiver(X,Y,P)

