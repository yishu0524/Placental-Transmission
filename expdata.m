%load the experimental data
plasma = readtable('plasma data.csv');
%%extract data
Monkey1=plasma(plasma.AnimalID==40690,:);
Monkey2=plasma(plasma.AnimalID==34795,:);
Monkey3=plasma(plasma.AnimalID==45010,:);
Monkey1=Monkey1(:, {'days', 'viralLoad'})
Monkey2=Monkey2(:, {'days', 'viralLoad'})
Monkey3=Monkey3(:, {'days', 'viralLoad'})
%%
Monkey1=Monkey1(~any(ismissing(Monkey1),2),:);
Monkey2=Monkey2(~any(ismissing(Monkey1),2),:);
Monkey3=Monkey3(~any(ismissing(Monkey1),2),:);
T_end =table(166,1,'VariableNames',{'days', 'viralLoad'})
Monkey1=[Monkey1;T_end]
Monkey2=[Monkey2;T_end]
%%
T=166;
l=2; %width of the placenta
D=0.0069;
nx=2000;
nt=1000;

x = linspace(0,l,nx); %2–2.5 cm (0.8–1 inch) in thickness, 
t = linspace(0,T,nt); %if we run simulation for 20 days
m = 0;
%%
Monkey1_VL=interp1(Monkey1.days,Monkey1.viralLoad,t,'pchip','extrap');
Monkey2_VL=interp1(Monkey2.days,Monkey2.viralLoad,t,'pchip','extrap');
Monkey3_VL=interp1(Monkey3.days,Monkey3.viralLoad,t,'pchip','extrap');
%%
figure(5)
semilogy(Monkey1.days,Monkey1.viralLoad,'o-')
hold on
semilogy(Monkey2.days,Monkey2.viralLoad,'o-')
hold on
semilogy(Monkey3.days,Monkey3.viralLoad,'o-')
hold off
title('Viral Load Data from mother Macaques')
xlabel('Time t (days)')
ylabel('CMV DNA Viral Load copies/ml')
legend('40690','34795','45010')
%%
%plot the 1d interpolation 

figure(1)
plot(t,Monkey2_VL) 
hold on
plot(Monkey2.days,Monkey2.viralLoad,'o')

hold off
title('Viral Load Data from mother Macaque 34795')
xlabel('Time t (days)')
ylabel('CMV DNA Viral Load copies/ml')
%legend('40690','34795','45010')
%%
s1=t;
s2=Monkey1_VL;
p=2/D;


options=odeset('RelTol',1e-4,'AbsTol',1e-4,'NormControl','off','InitialStep',1e-7) 
sol = pdepe(0,@unsatpde,@unsatic,@unsatbc,x,t,options,s1,s2,p,D); 
u= sol(:,:,1);
%%
placenta_growth=10000./(1+ exp(-(0.05*(t-80))));
%%
[ux,ut]=gradient(u,l/nx,T/nt);

virus_numeric=trapz(ux(:,1) .*transpose(placenta_growth) )*T/nt*D

%%
figure(2)
plot(t,placenta_growth)
title('Placental growth over pregancy')
xlabel('Time t (days)')
ylabel('Placental Surface Area (mm^2)')
legend('Placental Growth Function')
%%
Q = ux(:,1).* D.*transpose(placenta_growth);
virus_total=trapz(t,Q)
avg=trapz(t,placenta_growth)/T;
%%
figure (4)
plot(t,Q,'b')
%hold on
%plot(t, ux(:,1).* D*avg,'g')
%hold off
ylim([0 inf])
title('Viral flux on the infant of mother Macaque 40690 over time')
xlabel('Time t (days)')
ylabel('CMV DNA Viral Load copies/(μl*mm*day)')
%legend('Growth Function','Average Size')
%%
figure (3)
surf(x,t,u,'FaceAlpha',0.5,'EdgeColor','none')
hold on 
scatter3([2,2,2,2,2,2,2],Monkey2.days,Monkey2.viralLoad,'or','filled')
xlabel('Distance x (mm)')
ylabel('Time t (days)')
ylim([0,166])
zlabel('CMV DNA Viral Load copies/μl')
ax=gca;
ax.YDir="reverse"
hold off
title('Placental Tranport of CMV for mother Macaque 34795')

%%
t=transpose(t)
%%
function [c,f,s] = unsatpde(x,t,u,DuDx,s1,s2,p,D) 
c =1/D;% 1/(1.8*10^-3) ; %inverse of the diffusion coefficent
f = DuDx;
s =-p *u; %death of the cmv %death rate is 1, 2 virus day per day
end 
% ------------------------------------------------------------------------- 
function u0 = unsatic(x,s1,s2,p,D) 
u0 = 0; %initial is 0
end 
% ------------------------------------------------------------------------- 
function [pl,ql,pr,qr] = unsatbc(xl,ul,xr,ur,t,s1,s2,p,D) 
pl = ul; %u(l,t)=1
ql = 0; 
pr = ur-interp1(s1,s2,t,'linear','extrap'); %u(r,t)=0
qr = 0;
end 