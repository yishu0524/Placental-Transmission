for i = 2:27
    experimental_data('animal_data_cleared.csv', i)
end
%%
function experimental_data(filename, column_id)
%read data
    plasma = readtable(filename);

    %%
    %%select monkey and clean data
    monkey_data = plasma(:,[1,column_id]);
    monkey_data=monkey_data(~any(ismissing(monkey_data),2),:);
    monkey_ID=char(monkey_data.Properties.VariableNames(2));
    monkey_ID = strrep(monkey_ID,'_','-');
    monkey_ID = strcat(int2str(column_id-1),monkey_ID(1:end-5));
    monkey_data.Properties.VariableNames{1}='days';
    monkey_data.Properties.VariableNames{2}='viralLoad';
    %%
    %%parameter for pde
    T = table2array(monkey_data(end,1));
    if T < 100
        T_end =table(100,0,'VariableNames',{'days', 'viralLoad'});
        monkey_data=[monkey_data;T_end];
    end
    if T < 166
        T_end =table(166,0,'VariableNames',{'days', 'viralLoad'});
        monkey_data=[monkey_data;T_end];
        T = 166;
    end
        
    l=2; %width of the placenta
    D=0.01;
    nx=1000;
    nt=1000;

    x = linspace(0,l,nx); %2–2.5 cm (0.8–1 inch) in thickness, 
    t = linspace(0,T,nt); %if we run simulation for 20 days
    m = 0;
    %%
    monkey_data_VL=interp1(monkey_data.days,monkey_data.viralLoad,t,'pchip','extrap');
    %%
    figure1=figure(1)

    plot(monkey_data.days,monkey_data.viralLoad,'o')
    hold on
    plot(t,monkey_data_VL)
    ylim([0 inf])
    xlim([-inf 166])
    hold off
    title(monkey_ID)
    xlabel('Time t (days)')
    ylabel('CMV DNA Viral Load copies per microliter')
    legend('Experimental Data','pchip interpolation')
    saveas(figure1, strcat(monkey_ID,'.png'))
    %%
    %% pass the monkey data into pde
    %%
    s1=t;
    s2=monkey_data_VL;
    p=0.1/D;


    options=odeset('RelTol',1e-4,'AbsTol',1e-4,'NormControl','off','InitialStep',1e-7) 
    sol = pdepe(0,@unsatpde,@unsatic,@unsatbc,x,t,options,s1,s2,p,D); 
    u= sol(:,:,1);
    %%
    placenta_growth=10000./(1+ exp(-(0.05*(t-80))));
    %%
    [ux,ut]=gradient(u,l/nx,T/nt);

    virus_numeric=trapz(ux(:,1) .*transpose(placenta_growth) )*T/nt*D

    %%
    figure2=figure(2)
    plot(t,placenta_growth)
    ylim([0 11000])
    xlim([0 166])
    title('Placental growth over pregancy')
    xlabel('Time t (days)')
    ylabel('Placental Surface Area (mm^2)')
    legend('Placental Growth Function')
    saveas(figure2, 'placenta_growth.png')
    %%
    Q = ux(:,1).* D.*transpose(placenta_growth);
    virus_total=trapz(t,Q)
    avg=trapz(t,placenta_growth)/T;
    %%
    figure3=figure (3)
    plot(t,Q,'b')
    %hold on
    %plot(t, ux(:,1).* D*avg,'g')
    %hold off
    ylim([0 inf])
    xlim([0 166])
    title(monkey_ID)
    xlabel('Time t (days)')
    ylabel('CMV DNA Viral Load copies/(microliter*mm*day)')
    saveas(figure3, strcat(monkey_ID,'flux.png'))
    %%
    figure4= figure (4)
    surf(x,t,u,'FaceAlpha',0.5,'EdgeColor','none')
    hold on 
    scatter3(repelem(2,size(monkey_data,1)),monkey_data.days,monkey_data.viralLoad,'or','filled')
    xlabel('Distance x (mm)')
    ylabel('Time t (days)')
    ylim([0,166])
    zlim([0,inf])
    zlabel('CMV DNA Viral Load copies per microliter')
    ax=gca;
    ax.YDir="reverse"
    hold off
    title(monkey_ID)
    saveas(figure4, strcat(monkey_ID,'flow.png'))
    %%
    t=transpose(t);
    Table=table(t,Q);
    Table.Properties.VariableNames{1}='days';
    %%
    Table.Properties.VariableNames{2}=monkey_ID;
    writetable(Table, strcat(monkey_ID,'flow.csv'))
%%
end
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