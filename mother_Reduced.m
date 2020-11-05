
function dydt = mother_Reduced(t,y)
%Units: time(day) length(mm)

V=y(1); %Viral load
E=y(2); %virus-specific immune effector cells
RI=y(3); % actively-infected cells

%parameters
rho=5; % Virion induced immune response cells/(virions*day)
n=50; %Productivity of infected cell
delta=0.2; %Rate of viral induced cell death
c=0.3; %Rate of viral clearance
k=1e-4; %Infection rate constant
m=0.1; %Immune induced cell lysis
a0=0.2; %Exit and reactivation rate for monocytes
k=2e-3; %rate of latency
lambdaE=4e-2; %Homeostatic replenishment of immune cells
lambda=1e-3; %Cell replenishment rate
effector=9; %HCMV-specific effector cell term
epsilonS=0; %Level of immune suppression
f=1; %Number of infecting virions per cell
rS=4e2; %Equilibrium level of susceptible cells
rL=4e-2; %Equilibrium level of latent cells
tH=2; %Half-life of virions during antiviral treatment
tD=1.5; %Doubling time of the virions
Ee=10; %Equilibrium level of HCMV-specific effector cells
Ve=1e-2; %Equilibrium level of virions


Vderivative=n*delta*RI-c*V-f*k*rS*V;
Ederivative=(lambdaE*(1-E/effector)*E + rho *V);
RIderivative=k*rS*V-m*E*RI;

dydt = [Vderivative; Ederivative; RIderivative];

end