%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FINITE DIFFERENCE DISCRETIZATION OF THE TWO DIMENSIONAL NAVIER STOKE'S EQUATION 
%                                                        CAVITY FLOW
%                                                        PROBLEM (EXPLICIT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISCRETIZATION: : SECOND ORDER CENTRAL DIFFERENCE SCHEME  FOR DIFFUSION
%                             SECOND ORDER CENTRAL DIFFERENCE SCHEME  FOR  PRESSURE CORRECTION
%                             FIRST ORDER BACKWARD DIFFERENCING FOR CONVECTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE WRITTEN BY: NAHOM ALEMSEGED WORKU
% GENERAL EQUATION IN VECTOR FORM: DU/DT + (U . DIV) U = -1/RHO * GRAD(P) + DIV(GRAD(U))
% BOUNDARY CONDITION: U,V = 0 @ X=0,Xp, Y_BOTTOM, U=1 AND V = 0 FOR Y_TOP
% DP/DX = 0 @ X = 1,Xp, DP/DY = 0 @Y =BOTTOM,P = 0 @ Y_TOP                      
% DOMAIN = [0 X Xp] [0 X Yp]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLE DECLERATION
Xp = 8;    Yp = 8;      Tp = 10;
nx = 9;    ny = 9;        nt = 11;
dx = Xp/(nx-1);    dy = Yp/(ny - 1);        dt = Tp/(nt-1); 
P = zeros(nx,ny); 
b = zeros(nx,ny);
u = zeros(nx,ny);                   un = zeros(nx,ny);                    
v = zeros(nx,ny);                   vn = zeros(nx,ny);  pn = zeros(nx,ny); 
X = [0:dx:Xp];          Y = [0:dy:Yp];
maxit = 50;            rho = 1;        vis = 0.1;
% u(1,:) = 1 ;
% BOUNDARY CONDITIONS AND RIGHT HAND SIDE OF PRESSURE TERM
for k = 1:nt
    un = u;
    vn = v;
    for i=2:nx-1
        for j = 2:ny-1
            b(i,j) = (-rho.*(((un(i+1,j) - un(i-1,j)) / (2*dx))^2) + (((un(i+1,j) - un(i-1,j)) / (2*dx)) * ((vn(i,j+1) - vn(i,j+1)) / (2*dy))) +  ((vn(i,j+1) - vn(i,j+1)) / (2*dy)) + (rho * (((un(i+1,j) - un(i-1,j)) / (2*dx)) + ((vn(i,j+1) - vn(i,j+1)) / (2*dy)))/dt));                    
        end
    end
    %    FIND THE CORRECTED PRESSURE TERMS
    for it = 1:maxit
        Pd = P;
        for i=2:nx-1
            for j=2:ny-1
                P(i,j) = (((dy^2).*(P(i-1,j) + P(i+1,j))) +  ((dx^2).*(P(i,j-1) + P(i,j+1))) - (b(i,j) * (dx^2) * (dy^2))) / (2.*((dx^2) + (dy^2)))
            end
        end
    end
% COMPUTATION OF NAVIER STOKES EQUATION
        P(:,1) = P(:,2);        P(:,(end)) = P(:,end - 1);
        P(1,:)  = P(2,:);     P(end,:) = 0;
     
        un = u;     vn = v;
        for i = 2:nx - 1
            for j = 2:ny - 1
                u(i,j) = un(i,j) + dt*((-un(i,j) * ((un(i,j) - un(i-1,j))/dx)) - (vn(i,j) * ((un(i,j) - un(i,j-1))/dy)) - ((1/rho) * (P(i-1,j) - P(i-1,j))/(2*dx)) + (vis * ((un(i-1,j) - (2*un(i,j))+ un(i+1,j))/(dx^2))) + (vis * ((un(i,j-1) - (2*un(i,j))+ un(i,j+1))/(dy^2))));
                v(i,j) = vn(i,j) + dt*((-un(i,j) * ((vn(i,j) - vn(i-1,j))/dx)) - (vn(i,j) * ((vn(i,j) - vn(i,j-1))/dy)) - ((1/rho) * (P(i-1,j) - P(i-1,j))/(2*dx)) + (vis * ((vn(i-1,j) - (2*vn(i,j))+ vn(i+1,j))/(dx^2))) + (vis * ((vn(i,j-1) - (2*vn(i,j))+ vn(i,j+1))/(dy^2))));
            end
        end
        u(:,1) = 0;     u(:,end) = 0;       u(end,:) = 1;       u(1,:) = 0;
        v(:,1) = 0;     v(:,end) = 0;       v(end,:) = 0;       v(1,:) = 0;
end
% VISUALIZATION OF PRESSURE AND VELOCITY
figure(1)       
contourf(X,Y,P)
colorbar
xlabel('X')
ylabel('Y')
title('PRESSURE PROFILE')
grid on;
figure(2)
contourf(X,Y,u)
colorbar
xlabel('X')
ylabel('Y')
title('VELOCITY(U) PROFILE')
grid on;
% figure(3)
% contourf(X,Y,v)
% colorbar
% grid on;

    










