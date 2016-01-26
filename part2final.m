% FTCS diffusion model
% written by agt 1/22

clear all
figure(1)
clf

%% initialize - cape thompson

%load cape thompson data
load cape_thompson_copy.dat;
zdata = cape_thompson_copy(:,1); %depth column
Tdata = cape_thompson_copy(:,2); %temperature column

% find the equation of the line that fits the bottom straight part of cape
% thompson data
zlast = zdata(end); %define last data point
Tlast = Tdata(end);
slope = diff(Tdata)./diff(zdata); %find slope between each data point
dTdz_data = mean(slope(end-10:end)); %average of slopes between the last ten points
Told = Tlast -(dTdz_data*zlast); %solve equation of line for intercept at surface
T0_data = Told + (dTdz_data*zdata); %final equation of straight line

%% initialize - model

%constants
k = 2; %thermal cond w/mK
rho = 2000; % density kg/m3
Cp = 2000; % heat capacity 3/kg K

%set up depth array
zmax = 400; %maximum depth in m
dz = 5; %spacing in m
z = dz/2:dz:zmax-(dz/2); %depth array
N = length(z); %number of 'boxes', one surrounding each depth in z array

%set up time array
nyears = 100;
tmax = nyears*24*3600*365; %max time in seconds, nyears long
ndays = 20;
dt = ndays*24*3600; %increment of time step in s, ndays long
t = 0:dt:tmax; %creates time array

%setting up numerical solution
q = zeros(N,1); %creates heat flux array full of zeros
T = zeros(N,1); %creates a temp array full of zeros
dTdz0 = dTdz_data; %rename previously calculated slope of data
T = Told+(dTdz0*z); %rewrites temp array as line that fits straight part of data, starting point for solution
T0=T; %saves above array as T0 so that T can be rewritten
Ts_new = -5.2; %new surface temperature after temp step change
dTdz(N) = dTdz_data; % prescribed heat flux at the bottom in K/m, same as data

imax = length(t);
nplots=100; %number of plots generated
tplot = tmax/nplots; %time between each plot

%% run

for i = 1:imax
    
    T(1) = Ts_new; %change surface temp to new temp and hold it there
    
    dTdz(1:N-1) = diff(T)/dz; %calculate T gradient between each box
    q = -k*dTdz; %calculate heat flux
    dqdz = diff(q)/dz; % rate of change of temp T between each box

    % update T
    T(2:N) = T(2:N) - (1/(rho*Cp))*dqdz*dt; %updates all nodes except the top which we already set

    if(rem(t(i),tplot)==0) %only plot when (time in the loop)/(time between plots) doesn't have remainder
        figure(1)
        plot(T,z)
        hold on
        plot(T0,z,'r')
        plot(Tdata,zdata,'ko')
        
        xlabel('Temperature (°C)','fontname','arial','fontsize', 21)
        ylabel('Depth (m)', 'fontname', 'arial', 'fontsize', 21)
        set(gca, 'fontsize', 18, 'fontname', 'arial') 
        set(gca, 'YDIR', 'reverse') %change y-axis direction
        pause(0.1)
        hold off
    end

end

%% finalize

Tmod = interp1(z,T,zdata); %temperatures of model line at depth of each data point
chi2 = sum((Tdata - Tmod).^2); %chi squared test

if chi2<=.1 % test for goodness of fit, not actually sure what this number should be
    fprintf('good fit')
else
    fprintf('bad fit')
end


