%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FINITE DIFFERENCE  DISCRETIZATION OF TWO DIMENSIONAL BURGER FISHERS EQUATION 
%                                                         (UNSTEAD CONVECTION DIFFUSION EQUATION)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCHEME: IMPLICIT SHCEME
%DISCRETIZATION : FORWARD IN TIME , BACKWARD IN SPACE (UPWIND) FOR
%                   CONVECTION AND CENTRAL DIFFERENCE FOR DIFFUSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE  BY : NAHOM ALEMSEGED WORKU
% GENERAL EQUATION = dU/dt + U*dU/dX + V*dU/dY = VIS * (D2U/DX2 + D2U/DY2)
% AND 
% dV/dt + U*dV/dX + V*dV/dY = VIS * (D2V/DX2 + D2V/DY2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL CONDITIONS: U =2  0.5<=X<=1  AND  0.5<=Y<=1,
%                                 U = 1 EVERYWHERE ELSE
% BOUNDARY CONDITION ; U = 0, x = 0,2, AND Y = 0,Y = 2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DECLARE VARIABLES
Xp = 5; Yp = 5; ts = 4;
nx = 6;    ny = 6;    
nt =6;    dt = ts/(nt-1);    
dx = 5/(nx-1);   dy = 5/(ny-1);
vis = 1; 
u = zeros(nx,ny); %INITIALIZATION OF VELOCITY VECTOR
v = zeros(nx,ny);
X = [0:dx:Xp];  %DISRETIZATION IN X-DIRECTION
Y = [0:dy:Yp];  %DISCRETIZATION  IN Y-DIRECTION
t = [0:dt:ts];  %DISCRETIZATION OF TIME STEPS

% FILLING OUT MATRIX WITH INITIAL AND BOUDNARY CONDITIONS
for k=1:length(t)  
    if t(k) == 0 %LOOP TO HANDLE AND FILL INITIAL CONDITIONS
                 for i=1:(nx)
                      for j=1:(ny)
                            if ((X(i) >= 0 && X(i) < 1) || (Y(j) >= 0 && Y(j) < 1))
                                u(i,j) = 2;
                                v(i,j) = 1;
                                uv(i,j) = sqrt((u(i,j).^2) + (v(i,j).^2))
                           else
                                u(i,j) = 1;
                                v(i,j) = 1;
                                uv(i,j) = sqrt((u(i,j).^2) + (v(i,j).^2))
                          end
                      end
                 end
       else
        for i=1:nx   %LOOP TO HANDLE BOUNDARY CONDITIONS
            for j=1:ny
                 if ((X(i) == 0 || X(i) == X(end)) || (Y(j) == 0 || Y(j) == Y(end)))
                                u(i,j) = 0;
                                v(i,j) = 0;
                                uv(i,j) = sqrt((u(i,j).^2) + (v(i,j).^2))
            u    
            v
                 end
            end
        end
    end
end

 for it = 2:nt
        un = u; %ASSIGNING U VALUES FROM PREVIOUS TIME STEPS TO THE CURRENT
        vn = v;
        unvn(i,j) = sqrt((un(i,j).^2) + (vn(i,j).^2))
     for i=2:(nx-1)
         for j=2:(ny-1)
                   u(i,j) = ((((dt/dx) + (vis/(dx^2))) * (u(i-1,j))) + (((dt/dy) + (vis/(dy^2))) * u(i,j-1)) + un(i,j)) / (dt * ((u(i,j)/dx + v(i,j)/dy + (2*vis/(dx^2))) + (2*vis/(dy^2))))  
                   v(i,j) = ((((dt/dx) + (vis/(dx^2))) * (v(i-1,j))) + (((dt/dy) + (vis/(dy^2))) * v(i,j-1)) + vn(i,j)) / (dt * ((u(i,j)/dx + v(i,j)/dy + (2*vis/(dx^2))) + (2*vis/(dy^2)))) 
                   uv(i,j) = sqrt((u(i,j).^2) + (v(i,j).^2))
         end
     end
 end
%  VISULAIZATION OF TWO DIMENSIONAL PLOT
%  [X,Y] = meshgrid(X,Y)
   figure(1)
     contourf(X,Y,u,21,'LineStyle','none')
    colormap('jet')      
    colorbar
    hold off
     xlabel('X')
    ylabel('Y')
    title('U-VELOCITY COMPONENT')
     grid on;
        
  figure(2)
     contourf(X,Y,v,21,'LineStyle','none')
    colormap('jet')      
    colorbar
    hold off
     xlabel('X')
    ylabel('Y')
     grid on;
     title('V-VELOCITY COMPONENT')
     
      figure(3)
     contourf(X,Y,uv,21,'LineStyle','none')
    colormap('jet')      
    colorbar
    hold off
     xlabel('X')
    ylabel('Y')
     grid on;
     title('RESULTANT-VELOCITY COMPONENT')
      


 
     