%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SHOOTTING BISECTION METHOD FOR A TEMPERATURE PROFILE OF A RECTANGULAR FIN PROBLEM WITH
% INTERNAL HEAT GENERATION PARAMETER
%             CODE MODIFIED BY: NAHOM ALEMSEGED WORKU
%               TEMPERATURE PROFILE OF VARYING THERMAL CONDUCTIVITY
%%%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$%%

a=0; b=1;alpha=0; beta=1; tol=10^-6; m=1; n=1; Q=0.2;B= 0:1:2;G=0.2;
% aa=zeros(11,1);
  for k=1:3  
%     for j=1:11
        f=@(x,u) ([u(2);(((m.^2).*(u(1).^(n+1)))-((m.^2).*Q.*(1+(G.*u(1))))-(B(k).*(u(2).^2)))./(1+(u(1).*B(k)))]);
        %   f=@(x,u) ([(u(2)*(1+B(k).*u(1)));((m(j).^2)*(u(1).^(n+1))-((m(j).^2).*Q.*(1+(G.*u(1))))-(B(k).*(u(2)^2)))]);
        xL=-2;
        xU=2;%  Lower and upper values of the slope(Left and Right)
        imax=150;
            for i=1:imax
            xr=xL+0.5*(xU-xL);
            ic=[xr alpha];
            [x,u]=ode45(f,[0,1],ic);
            u(i)=u(end,1);
            err=u(i)-beta;
            if abs(err)<tol
                break;
            end
           if err>0
                xU=xr;  
                else
                    xL=xr;
           end
           end
                [x,u]=ode45(f,[0 1],ic);
                aa =-u(end,2);
             
           if k==1
            plot(x,u(:,1),'ok');  
            elseif k==2
                 plot(x,u(:,1),'*');
            elseif k==3
                  plot(x,u(:,1),'.');

           end  
                  hold on
           grid on
      end
   
xlabel('Thermo-Geometric Parameter (m)');
ylabel('(Nusselt No.)-d\theta /dt');
legend('\beta=0','\beta=1','\beta=2');
title('Plot  Fin Problem for n=1,\gamma=0.2 and Q=0.2')