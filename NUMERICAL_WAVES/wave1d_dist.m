%% Numerical Project: Waves in MITgcm
% This code solves the 2nd order 1D wave equation

% Run Parameters
L = 200*pi; % Half-domain length
dur = L*4; % Run time
L0 = L/10; % Decay scale of gaussian envelope
k0 = 20*pi/L0; % Wavenumber
dx = 2*pi/(20*k0); % Spatial step
c = 1; % Phase speed 
dt = dx/c; % Timestep
sig2 = (c*dt/dx)^2; % Non-dim. coefficient
w0=c*k0; % Frequency

% Initialize vectors
x = -L:dx:L; % Space
K = length(x);
t = 0:dt:dur; % Time
T = length(t);

% Set up matrices
Mdir = diag(-2*ones(1,K),0) + diag(ones(1,K-1),1) + diag(ones(1,K-1),-1);
Mdir(1,:)=0;
Mdir(end,:)=0;

ic=3;
switch ic
    case 1 %Carrier waves with gaussian envelope
        % Initial conditions
        phi = NaN(length(x),length(t));
        phi0 = exp(-x.^2/L0^2).*cos(k0*x);
        phid0 = (2*c*L0^(-2)*x.*cos(k0*x)+w0*sin(k0*x)).*exp(-x.^2/L0^2);
        phi(:,1)=phi0';
        phi(:,2)=phi(:,1)+dt*phid0'+1/2*sig2*Mdir*phi(:,1);

        for n=2:length(t)
            phi(:,n+1) = (2*eye(K)+sig2*Mdir)*phi(:,n)-phi(:,n-1);
        end

    case 2 % Gaussian pulse
        % Initial conditions
        phi = NaN(length(x),length(t));
        phi0 = exp(-x.^2/L0^2);
        phid0 = 2*c/L0^2*x.*phi0;
        phi(:,1)=phi0';
        phi(:,2)=phi(:,1)+dt*phid0'+1/2*sig2*Mdir*phi(:,1);

        for n=2:length(t)
            phi(:,n+1) = (2*eye(K)+sig2*Mdir)*phi(:,n)-phi(:,n-1);
        end
        
    case 3 % Gaussian pulse with neumann boundary conditions
        L = 200*pi;
        dur = L*8;
        L0 = L/10;
        k0 = 20*pi/L0;
        dx = L0/10;
        c = 1;
        dt = dx/c;
        sig2 = (c*dt/dx)^2;
        w0=c*k0;

        % Initialize vectors
        x = -L:dx:L;
        K = length(x);
        t = 0:dt:dur;
        T = length(t);

        % Set up matrices
        Mneu = diag(-2*ones(1,K),0) + diag(ones(1,K-1),1) + diag(ones(1,K-1),-1);
        Mneu(1,:)=0;
        Mneu(end,end-1)=2;

        % Initial conditions
        phi = NaN(length(x),length(t));
        phi0 = exp(-x.^2/L0^2);
        phid0 = 2*c/L0^2*x.*phi0;
        phi(:,1)=phi0';
        phi(:,2)=phi(:,1)+dt*phid0'+1/2*sig2*Mneu*phi(:,1);

        % Evaluate
        for n=2:length(t)
            phi(:,n+1) = (2*eye(K)+sig2*Mneu)*phi(:,n)-phi(:,n-1);
        end
end
%% Plot snapshots in time

for i=round(linspace(1,T,20))
close(gcf)
figure
    plot(x,phi(:,i))
    xlabel('x');
    ylabel('\phi');
    title('\phi as a function of x at different times');
    ylim([-2.01,2.01])
    set(gcf, 'Position',  [576, 252, 768, 576]) % presentation size
    saveas(gcf,'phi.png')
end


%%
i=1
figure
    plot(x,phi(:,i))
    xlabel('x');
    ylabel('\phi');
    tlabel = sprintf('t = %i',t(i));
    title(['\phi as a function of x at ' tlabel]);
    ylim([-2.01,2.01])
    xlim([-L/4 L/4]);
    set(gcf, 'Position',  [576, 252, 768, 576]) % presentation size
    saveas(gcf,'phi_packet_initial.png')

%%
tsel = [1 L/2/dt+4 L/dt+4 3*L/dt+4];
figure
    hold on;
    for i=tsel
        plot(x,phi(:,i))
    end
    xlabel('x');
    ylabel('\phi');
    legend('t = 0','t = 0.5 L','t = L','t = 3 L');
    title('\phi as a function of x at selected times');
    ylim([-2.01,2.01])
    xlim([x(1) x(end)]);
    set(gcf, 'Position',  [576, 252, 768, 576]) % presentation size
    saveas(gcf,'phi_packet.png')


%%
tsel = int16([1 L/2/dt+1 L/dt+1 5*L/2/dt+1 3*L/dt+1 7/2*L/dt+1]);
figure
    hold on;
    for i=tsel
        plot(x,phi(:,i))
    end
    xlabel('x');
    ylabel('\phi');
    legend('t = 0','t = 0.5 L','t = L','t = 2.5 L','t = 3 L','t = 3.5 L');
    title('\phi as a function of x at selected times');
    ylim([-2.01,2.01])
    xlim([x(1) x(end)]);
    set(gcf, 'Position',  [576, 252, 768, 576]) % presentation size
    saveas(gcf,'phi_pulse.png')


%% Make a movie
fps = 60;
mlength = 10;
vidfile = VideoWriter('phi_pulse.mp4','MPEG-4');
vidfile.FrameRate=fps;
open(vidfile);  
figure(1);
    set(gcf, 'Position',  [576, 252, 768, 576]) % presentation size
    frames = fps*mlength;
    for i=round(linspace(1,T,frames))
        plot(x,phi(:,i))
        xlabel('x');
        ylabel('\phi');
        tlabel = sprintf('t = %.1f L',t(i)/L);
        title(['\phi as a function of x at ' tlabel]);
        ylim([-2.01,2.01])
        xlim([x(1) x(end)]);
        drawnow;
        writeVideo(vidfile, getframe(gcf));
        disp(i);
     end
close(vidfile);
close(gcf)
