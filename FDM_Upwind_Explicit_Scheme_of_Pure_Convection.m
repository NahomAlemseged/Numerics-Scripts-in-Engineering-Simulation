%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FINITE DIFFERENCE  DISCRETIZATION OF LINEAR CONVECTION EQUATION WITH CONSTANT VELOCITY 
% SCHEME: IMPLICIT SHCEME
%DISCRETIZATION : FORWARD IN TIME AND BACKWARD IN SPACE
% AN UPPWIND SCHEME OF CONVECTION IS USED
% CODE  BY : NAHOM ALEMSEGED WORKU
% GENERAL EQUATION = dU/dt + C*dU/dX + C*dU/dY = 0
% INITIAL CONDITIONS: U =2  0.5<=X<=1  AND  0.5<=Y<=1,

%                                 U = 1 EVERYWHERE ELSE
% BOUNDARY CONDITION ; U = 0, x = 0,2, AND Y = 0,Y = 2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DECLARE VARIABLES
Xp = 5; Yp = 5; ts = 2;
nx = 6;    ny = 6;    
nt =6;    dt = ts/(nt-1);    
dx = 5/(nx-1);   dy = 5/(ny-1);
c = 1; 
u = zeros(nx,ny); %INITIALIZATION OF VELOCITY VECTOR
X = [0:dx:Xp];  %DISRETIZATION IN X-DIRECTION
Y = [0:dy:Yp];  %DISCRETIZATION  IN Y-DIRECTION
t = [0:dt:ts];  %DISCRETIZATION OF TIME STEPS
% FILLING OUT MATRIX WITH INITIAL AND BOUDNARY CONDITIONS
for k=1:length(t)  
    if t(k) == 0 %LOOP TO HANDLE AND FILL INITIAL CONDITIONS
                 for i=1:(nx)
                      for j=1:(ny)
                            if ((X(i) >= 0 && X(i) < 1) || (Y(j) >= 0 && Y(j) < 1))
                                u(i,j) = 2
                           else
                                u(i,j) = 1
                          end
                      end
                 end
   
    else
        for i=1:nx   %LOOP TO HANDLE BOUNDARY CONDITIONS
            for j=1:ny
                 if ((X(i) == 0 || X(i) == X(end)) || (Y(j) == 0 || Y(j) == Y(end)))
                                u(i,j) = 0;
            u    
                 end
            end
        end
    end
end

 for it = 2:nt
        un = u %ASSIGNING U VALUES FROM PREVIOUS TIME STEPS TO THE CURRENT
     for i=2:(nx-1)
         for j=2:(ny-1)
%              u(i,j) = un(i,j) - ((c*dt/dx)*(u(i,j) - u(i-1,j))) - ((c*dt/dy)*(un(i,j) - un(i,j-1))); %FILL EXPLICIT DISCRETIZED VALUES
                 u(i,j) = (un(i,j) + ((c*dt/dx)*(u(i-1,j))) + ((c*dt/dy)*(u(i,j-1))))/(1 + (c*dt/dx) + (c*dt/dy)) ; 
         end
     end
 end
%  VISULAIZATION OF TWO DIMENSIONAL PLOT
 [X,Y] = meshgrid(X,Y)
   figure(1)
     contourf(X,Y,u,21,'LineStyle','none')
     hold on;
%     contour(u)
    colormap('jet')      
    colorbar
    hold off
     xlabel('X')
    ylabel('Y')
     grid on;

 
     