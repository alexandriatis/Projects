%% Numerical Project: Waves in MITgcm
% This code solves the 1-D shallow water equations with linearized
% free surface MITgcm-style. 

% Run Parameters
L = 200*pi; % Half-domain length
L0 = L/10; % Decay scale for gaussian I.C.
g = 9.81; % Gravity
hmax = 100; % Max water depth
hmin = 5; % Min water depth
cg = sqrt(g*hmax); % Shallow water wave speed
dur = 10*L/cg; % Run time
dx = L0/10; % Spatial step
dt = 2.6*dx/cg; % Time step
dt = min(dt,dur/600); % Minimum time step for making a movie
epsg = 0.1; % Small stability parameter


% Initialize vectors
x = -L:dx:L; % Length dimension
K = length(x);
t = 0:dt:dur; % Time dimension
T = length(t);
h = hmax*ones(size(x)); % Bottom depth
for i=(K+1)/2:K
    h(i) = (hmin-hmax)/L*x(i)+hmax;
end

% Initial conditions
eta = zeros(K,T); % Sea surface height
eta0 = exp(-(x+L/4).^2/L0^2); % Gaussian initial condition
eta(1:(K+1)/2,1) = eta0(1:(K+1)/2)';
u = zeros(K,T); % Water speed
u0 = cg*eta0./h; % Gaussian initial condition
u(1:(K+1)/2,1) = u0(1:(K+1)/2)';
eta(1,1)=0; % removes small leftover from exponential 
u(1,1)=0; % removes small leftover from exponential 

% Define matrices

H1 = zeros(K);
for k=2:K-1
    H1(k,k) = h(k+1)-h(k-1);
end
H1(1,1)=2*(h(2)-h(1));
H1(end,end)=2*(h(end)-h(end-1));

H2 = zeros(K);
for k=2:K-1
    H2(k,k-1) = -h(k);
    H2(k,k+1) = h(k);
end
H2(1,1) = -2*h(1);
H2(1,2) = 2*h(1);
H2(end,end-1) = -2*h(K);
H2(end,end) = 2*h(K);

M2h = zeros(K);
for k=2:K-1
    M2h(k,k-1) = h(k);
    M2h(k,k) = -2*h(k);
    M2h(k,k+1) = h(k);
end
M2h(1,1)=h(1);
M2h(1,2)=-2*h(1);
M2h(1,3)=h(1);
M2h(end,end-2)=h(end);
M2h(end,end-1)=-2*h(end);
M2h(end,end)=h(end);

M1h = zeros(K);
for k=2:K-1
    M1h(k,k-1) = h(k-1)-h(k+1);
    M1h(k,k+1) = h(k+1)-h(k-1);
end
M1h(1,1)=4*(h(1)-h(2));
M1h(1,2)=4*(h(2)-h(1));
M1h(end,end-1)=4*(h(end-1)-h(end));
M1h(end,end)=4*(h(end)-h(end-1));

M1 = diag(1*ones(1,K-1),1) + diag(-1*ones(1,K-1),-1);
M1(1,1)=0;
M1(1,2)=0;
M1(end,end-1)=0;
M1(end,end)=0;

% Evaluation

G = zeros(K,2);
for n=1:T-1
    tildeh = (h+eta(:,n));
    % Equation 69
    Gadv = zeros(K,1);
    for k=2:K-1
        Gadv(k) = 1/(2*dx*tildeh(k))*u(k,n)*(tildeh(k+1)*u(k+1,n)-tildeh(k-1)*u(k-1,n));
    end

    G(:,2) = -Gadv;
    if n==1
        G12 = G(:,2);
    else
        G12 = (3/2+epsg)*G(:,2)-(1/2+epsg)*G(:,1);
    end
    G(:,1)=G(:,2);
    ustar = u(:,n)+dt*G12;
    
    % Equation 72
    etastar = eta(:,n) - dt/(2*dx)*(H1*ustar+H2*ustar);
    
    % Equation 67
    eta(:,n+1) = -1/dt^2*(g/dx^2*M2h + g/(4*dx^2)*M1h-1/dt^2*eye(K))^(-1)*etastar;
    
    % Equation 75
    u(:,n+1) = ustar - g*dt/(2*dx)*M1*eta(:,n+1);
end


% Plot final state
i=n;
figure
    plot(x,eta(:,i));
    hold on
    plot(x,u(:,i));
    plot(x,etastar);
    plot(x,ustar);
    plot(x,eta(:,i+1));
    plot(x,u(:,i+1));
    plot(x,-2*h/hmax,'k--');
    xlabel('x');
    ylabel('\eta');
    tlabel = sprintf('t = %.1f L/cg',t(i)/(L/cg));
    title(['\eta as a function of x at ' tlabel]);
    ylim([-2.6,2.6])
    legend('\eta','u','\eta *','u *','\eta+1','u+1','2h/h_{max}');
    set(gcf, 'Position',  [576, 252, 768, 576]) % presentation size
    saveas(gcf,'eta.png')

%% Plot snapshots in time
tsel = [1 L/(cg*dt) 2*L/(cg*dt) 3*L/(cg*dt) 4*L/(cg*dt)];
tsel = round(tsel);
figure
    hold on;
    for i=tsel
        plot(x,eta(:,i))
    end
    plot(x,-2*h/hmax,'k--');
    xlabel('x');
    ylabel('\eta');
    legend('t = 0','t = L/c_{g}','t = 2L/c_{g}','t = 3 L/c_{g}','t = 4 L/c_{g}','2h/h_{max}');
    title('\eta as a function of x at selected times');
    ylim([-2.6,2.6])
    xlim([x(1) x(end)]);
    set(gcf, 'Position',  [576, 252, 768, 576]) % presentation size
    saveas(gcf,'eta_pulse.png')

%% Make a movie
fps = 60;
mlength = 10;
vidfile = VideoWriter('eta_pulse.mp4','MPEG-4');
vidfile.FrameRate=fps;
open(vidfile);  
figure(1);
    set(gcf, 'Position',  [576, 252, 768, 576]) % presentation size
    frames = fps*mlength;
    for i=round(linspace(1,T,frames))
        clf;
        plot(x,eta(:,i))
        hold on;
        plot(x,-2*h/hmax,'k--');
        xlabel('x');
        ylabel('\eta');
        tlabel = sprintf('t = %.1f L/cg',t(i)/(L/cg));
        title(['\eta as a function of x at ' tlabel]);
        legend('\eta','2h/h_{max}');
        ylim([-2.6,2.6])
        xlim([x(1) x(end)]);
        drawnow;
        writeVideo(vidfile, getframe(gcf));
        disp(i);
     end
close(vidfile);
close(gcf)

%% Check conservation rules
mass = NaN(1,T);
momentum = NaN(1,T);
ke = NaN(1,T);
pe = NaN(1,T);
energy = NaN(1,T);
for n=1:T
    tildeh = (h(:)+eta(:,n));
    mass(n) = sum(tildeh);
    momentum(n) = sum(tildeh.*u(:,n));
    ke(n) = sum(1/2*tildeh.*(u(:,n).^2));
    pe(n) = sum(1/2*g*eta(:,n).^2);
    energy(n) = ke(n)+pe(n);
end

figure
    plot(t,mass/mass(1));
    hold on
    plot(t,momentum/momentum(1));
    plot(t,ke/ke(1));
    plot(t,pe/pe(1));
    plot(t,energy/ke(1));

    xlabel('t');
    xticks([1 2*L/(cg) 4*L/(cg) 6*L/(cg) 8*L/(cg)]);
    xticklabels({'0','2L/c_{g}', '4L/c_{g}', '6L/c_{g}','8L/c_{g}'})
    ylabel('conserved quantity');
    title('Conserved quantities as a function of time');
    xlim([t(1),t(end)]);
    legend('mass (h+\eta)','momentum ((h+\eta)u)','KE (1/2(h+\eta)u^2)','PE (1/2g(h+\eta)^2)','E = PE+KE');
    set(gcf, 'Position',  [576, 252, 768, 576]) % presentation size
    saveas(gcf,'conservation.png')    
