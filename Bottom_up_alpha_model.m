    % Matlab code for implementing and plotting the model presented in Karvat 
% and Landau (2023), "A role for bottom-up alpha oscillations in temporal 
% integration", Journal of Cognitive Neuroscience.
%
% Version 1.0, Golan Karvat, 2023

clearvars;
figure(627); clf; hold on
siz = get(0,'screensize');
pos = [siz(3)*3/10, siz(4)/10, siz(3)*4/10, siz(4)*8/10];
set(gcf,'position',pos); 
rows = 4;
cols = 2;

% Parameters and constants
thr = exp(-1);      % threshold underwhich an input generates a seperate percept
ISI = 75;           % Inter stimuuls interval, in ms
dt = 0.001;         % time resolution of the plot
ISI = ISI*dt;       % ISI translated into time resolution units.
t = 0:dt:0.5;       % time vector
TAU = [1,2]*30*dt;  % time decay constant 
t0 = 100*dt;        % time in which the first "flash" appears (100 ms).
A = [1, 1];         % Amplitude of the impulse response

% Configuration structure for the oscillation-pulse interaction (panels C-H);
CFG         = [];
CFG.thr     = 0;    % below this threshold, the response effect is damped
CFG.perThr  = exp(-1); % area to fill as "effective"
CFG.Fosc    = 10;   % frequency of the oscillation
CFG.Aosc    = 1;    % amplitude of the oscillation (rel. impulse response)
CFG.Aosc2   = 0;    % amplitude of the burst oscillation oscillation (rel. impulse response)
CFG.Tosc2   = 0;    % Timing of the burst (rel impulse response).
CFG.phi     = 0*pi/4; % phase shift in radians
CFG.ISI     = 50;   % inter-stimuls intereval, in ms

tits = {'A','B','C','D','E','F','G','H'};   % Title strings

for ti = 1:length(TAU) % loop through time constants, resembling persistence and affected by lateral inhibition
    subplot(rows,cols,ti); cla; hold on
    tau = TAU(ti);
    RSP = exp(1-tau*t/dt).*t/tau; % The response, modeled as a dampened impulse response following the "ERP" in van Diepen & Mazaheri, 2018
    RSP1 = [zeros(1,t0/dt),RSP(1:end-t0/dt)]/max(RSP)*A(ti); % Response to the first flash
    RSP2 = [zeros(1,round((t0+ISI)/dt)),RSP(1:end-round((t0+ISI)/dt))]/max(RSP)*A(ti); % Response to the second flash, first ISI
    RSP3 = [zeros(1,round((t0+2*ISI)/dt)),RSP(1:end-round((t0+2*ISI)/dt))]/max(RSP)*A(ti); % Response to the second flash, second ISI
    
    [m1,i1] = max(RSP1);
    [m2,i2] = max(RSP2);
    [m3,i3] = max(RSP3);
    
    % plot the first response
    curr = RSP1;    
    pcurr = curr(curr>=thr);
    pcurr([1 end]) = thr;
    patch(t(curr>=thr),pcurr,ones(1,3)*0.65,'LineStyle','none');
    plot(t,curr,'k','linewidth',2)
    scatter(t(i1),m1,100,'b','fill');
    
    % plot the second response with first ISI, and mark it if "detected" 
    % (= arrived when the first response is already below threshold)
    curr = RSP1+RSP2;
    curr = curr/max(curr)*A(ti);
    impPt = round((t0+ISI)/dt);
    curr(1:impPt+3) = RSP1(1:impPt+3);
    plot(t,curr,'--k','linewidth',2);    
    if ti==2
        scatter(t(i2-3),m2,100,'r','fill');
    end
    
    % plot the second response with second ISI
    curr = RSP1+RSP3;
    curr = curr/max(curr)*A(ti);
    impPt = round((t0+2*ISI)/dt);
    curr(1:impPt+1) = RSP1(1:impPt+1);
    plot(t,curr,':k','linewidth',2)    
    scatter(t(i3),m3,100,'r','fill');
    
    tt(ti) = title(tits{ti});
    tt(ti).Units = 'normalized';
    tt(ti).Position(1) = 0; 
    ylim([0 1])
    axis off
end

% panels C-H
for spi = 1:6 
    cfg{spi} = CFG;
    cfg{spi}.tit = tits{spi+2};    
end
% Define varying parameters in each panel
cfg{2}.Aosc     = 0.25;
cfg{3}.phi      = -2*pi/4;
cfg{4}.phi      = -2*pi/4;
cfg{4}.Fosc     = 12;
cfg{5}.ISI      = 20;
cfg{6}.ISI      = 20;
cfg{6}.Aosc     = 0.25;
cfg{6}.Aosc2    = 1.5;
cfg{6}.Tosc2    = 10; 

for spi = 1:6
    subplot(rows,cols,spi+2); hold on
    OscRsp(cfg{spi}); % function that plots the osciilation and impulse response,
                      % and determines if the second pulse is detected.
end

set(findall(gcf,'-property','fontsize'),'fontsize',18,'fontweight','bold')

function OscRsp(cfg)
% get the values from the configuration structure
thr = cfg.thr;
Aosc = cfg.Aosc;
Fosc = cfg.Fosc;
phi = cfg.phi;
ISI = cfg.ISI;
perThr = cfg.perThr;

% constants
dt = 0.001;
ISI = ISI*dt;
t = 0:dt:0.5;
tau = 30*dt;
t0 = 100*dt;
A = 1;

% Model the impulse response (similar to above)
RSP = exp(1-tau*t/dt).*t/tau;
RSP1 = [zeros(1,t0/dt),RSP(1:end-t0/dt)]/max(RSP)*A;
RSP2 = [zeros(1,round((t0+ISI)/dt)),RSP(1:end-round((t0+ISI)/dt))]/max(RSP)*A;
Osc = Aosc*sin(2*pi*Fosc.*(t-t0)+phi); % the (alpha) oscillation
cyc = round(1/(Fosc*dt)); % one cycle duration
if cfg.Aosc2 % generate a "burst" with higher freuqncy
    Fosc2 = Fosc*4;
    Osc2 = cfg.Aosc2*sin(2*pi*Fosc2.*(t-t0)+phi);
    tb = round((t0+tau+ISI)/dt)+cfg.Tosc2; % time of the center of the burst
    gauss = zeros(size(t));
    bToPlot = nan(size(t));
    Hg = 50; % half the gaussian
    gauss(tb+[-Hg:Hg]) = gausswin(Hg*2+1);    
    burst = Osc2.*gauss;     
    Osc = Osc+burst;
    bToPlot(tb+[-Hg:Hg]) = burst(tb+[-Hg:Hg]);
    cyc = round(1/(Fosc2*dt));
end

curr = RSP1+Osc;
t1 = find(curr>Aosc,1,'first'); % the first point in which the impulse has an excitatory effect above the oscillation amplitude
dumpt = find (t>(t(t1)) & curr<=thr,1,'first'); % the first timepoint in whic hsignal went below a shreshold
RSP1(dumpt:end)=0; % set the response from this point onward to 0 
% (based on population refractory periods, (Buzs?ki & V?r?slakos, 2023; Sanchez-Vives & McCormick, 2000))
ind(1) = find(t>=t0 & curr>=perThr,1,'first'); % the first point of the "effective" (persistence) area

% repeat with the response to the second "flash"
curr = RSP1+RSP2+Osc;
t2 = find(RSP2+Osc>Aosc,1,'first');
dumpt = find (t>t(t2) & curr<=thr,1,'first');
RSP2(dumpt:end)=0;
curr = RSP1+RSP2+Osc;
ind(2) = find(t>t(ind(1)) & curr<=perThr,1,'first'); % the last point of the "effective" (persistence) area
inds = ind(1):ind(2);
pcurr = curr(inds); pcurr([1 end]) = perThr;
patch(t(inds),pcurr,ones(1,3)*0.65,'LineStyle','none'); % mark the "effective" area
try plot(t, bToPlot-0.7,'linewidth',3,'color',ones(1,3)*0.5); end % plot the burst, if exists
plot(t,curr,'k','linewidth',2)

% mark the peak of the response to the first flash
[~,I1] = max(RSP1);
scatter(t(I1), curr(I1),75,'b','fill');

% mark the peak of the response to the first flash, only if detected
[m2,I2] = max(curr(ind(2):ind(2)+cyc));
I2 = I2+ind(2)-1;
if m2>max([perThr,Aosc]); scatter(t(I2),curr(I2),75,'r','fill'); end
ys = get(gca,'ylim');
% mark the time of arriveal of te sceonds response
plot([0 0]+t0+ISI+tau/2,ys,'--','color',ones(1,3)*0.,'linewidth',2);

% title
tt = title(cfg.tit);
tt.Units = 'Normalized'; tt.HorizontalAlignment = 'left';
tt.Position(1) = 0; tt.FontSize = 16;
ylim([-2 2.75]);
axis off
end