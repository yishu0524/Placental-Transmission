%viral load on mother side
initial_values = [1e-4, 0, 0,4e2,0]; % (V, E, RI, RS, RL) 
 [t_m,y] = ode15s(@mother,[0 280],initial_values); 
% xlabel('Time t');
% ylabel('Solution y');
% set(gca, 'YScale', 'log')
% legend('V')
% x = 0:.25:5;
% y = x.^2;
y_m=y(:,1)*1e3;
V = @(z) interp1(t_m,y(:,1)*1e3,z);  % Lookup table function

t2 = 0:1:280;  % Just to compare func to data.
figure (1)
plot(t_m,y(:,1)*1e3,'sb',t2,V(t2),'*r')

%%
l=2; %width of the placenta
T=280;
D=0.0069;
c2=mean(y_m);
nx=2000;
nt=1000;
x = linspace(0,l,nx); %2–2.5 cm (0.8–1 inch) in thickness, 
t = linspace(0,T,nt); %if we run simulation for 20 days
m = 0;
s1=t_m;s2=y_m;p=0.5/D;


options=odeset('RelTol',1e-4,'AbsTol',1e-4,'NormControl','off','InitialStep',1e-7) 
sol = pdepe(0,@unsatpde,@unsatic,@unsatbc,x,t,options,s1,s2,p,D); 
u= sol(:,:,1);
figure (2)
plot(x,u(end,:)) 
[ux,ut]=gradient(u,l/nx,T/nt);
virus_numeric=trapz(ux(:,end))*T/nt*D*28858
%%
figure (3)
surf(x,t,u)
title('Numerical solution')
xlabel('Distance x')
ylabel('Time t')
%%

function [c,f,s] = unsatpde(x,t,u,DuDx,s1,s2,p,D) 
c =1/D;% 1/(1.8*10^-3) ; %inverse of the diffusion coefficent
f = DuDx;
s =-p*u; %death of the cmv %death rate is 1, 2 virus day per day
end 
% ------------------------------------------------------------------------- 
function u0 = unsatic(x,s1,s2,p,D) 
u0 = 0; %initial is 0
end 
% ------------------------------------------------------------------------- 
function [pl,ql,pr,qr] = unsatbc(xl,ul,xr,ur,t,s1,s2,p,D) 
pl = ul-interp1(s1,s2,t,'linear','extrap'); %u(0,t)=1
ql = 0; 
pr = ur; %u(L,t)=0
qr = 0;
end 
%%