set(0, 'defaultaxesfontsize', 20);
set(0, 'DefaultFigureWindowStyle', 'docked');
set(0, 'DefaultLineLineWidth', 2);
set(0, 'Defaultaxeslinewidth', 2);

% Constants
c_c = 299792458; 
c_eps_0 = 8.8542149e-12; 
c_eps_0_cm = c_eps_0 / 100; 
c_mu_0 = 1 / c_eps_0 / c_c^2; 
c_q = 1.60217653e-19; 
c_hb = 1.05457266913e-34; 
c_h = c_hb * 2 * pi; 


RL = 0.9i; % Reflectivity of the left mirror
RR = 0.9i; % Reflectivity of the right mirror 

% Input pulse parameters
InputParasL.E0 = 1e5; % Amplitude of the electric field
InputParasL.we = 0;   % Initial energy (set to zero for this case)
InputParasL.t0 = 2e-12; % Initial time delay (s)
InputParasL.wg = 5e-13; % Pulse width (s)
InputParasL.phi = 0;   % Initial phase
InputParasR = 0;       % No right input initially

n_g = 3.5; % Group index of the waveguide
vg = c_c / n_g * 1e2; % Group velocity (cm/s)
Lambda = 1550e-9; % Wavelength of light (m)

plotN = 10; 

% Spatial and temporal grid
L = 1000e-6 * 1e2; % Waveguide length in cm
XL = [0, L]; % Spatial range for plots
YL = [0, InputParasL.E0]; % Range for electric field plots
Nz = 500; % Number of spatial points
dz = L / (Nz - 1); % Spatial step size (cm)
dt = dz / vg; % Time step size (s)
fsync = dt * vg / dz; % Scaling factor for field propagation

% Total number of time steps
Nt = floor(2 * Nz);
tmax = Nt * dt; % Maximum simulation time (s)
t_L = dt * Nz; % Time to travel the waveguide length (s)

% Spatial and time arrays
z = linspace(0, L, Nz).'; % Spatial grid points
time = nan(1, Nt); % Time array
InputL = nan(1, Nt); % Left input array
InputR = nan(1, Nt); % Right input array
OutputL = nan(1, Nt); % Left output array
OutputR = nan(1, Nt); % Right output array

% Electric fields
Ef = zeros(size(z)); % Forward-propagating field
Er = zeros(size(z)); % Backward-propagating field

% Input functions
Ef1 = @SourceFct; % Source function for forward field
ErN = @SourceFct; % Source function for reverse field

% Initialize time and inputs
t = 0;
time(1) = t;
InputL(1) = Ef1(t, InputParasL); % Initial left input
InputR(1) = ErN(t, InputParasR); % Initial right input

% Initialize outputs
OutputR(1) = Ef(Nz); % Right output (initially zero)
OutputL(1) = Er(1); % Left output (initially zero)

% Apply initial boundary conditions
Ef(1) = InputL(1); % Forward field at z=0
Er(Nz) = InputR(1); % Backward field at z=L

% Plot the initial fields
figure('name', 'Fields')
subplot(3, 1, 1)
plot(z * 10000, real(Ef), 'r'); % Plot forward field (real part)
hold off
xlabel('z(\mum)')
ylabel('E_f')

subplot(3, 1, 2)
plot(z * 10000, real(Er), 'b'); % Plot backward field (real part)
xlabel('z(\mum)')
ylabel('E_r')
hold off

subplot(3, 1, 3)
plot(time * 1e12, real(InputL), 'r'); hold on % Plot left input
plot(time * 1e12, real(OutputR), 'r--'); % Plot right output
plot(time * 1e12, real(InputR), 'b'); % Plot right input
plot(time * 1e12, real(OutputL), 'b--'); % Plot left output
xlabel('time(ps)')
ylabel('E')
hold off

% Main simulation loop
for i = 2:Nt
    t = dt * (i - 1); % Iteration of time
    time(i) = t;

    % Update inputs
    InputL(i) = Ef1(t, InputParasL); % Left input at time t
    InputR(i) = ErN(t, 0); % Right input at time t

    % Boundary conditions with mirrors added
    Ef(1) = InputL(i) + RL * Er(1); % Forward field at z=0 with left mirror
    Er(Nz) = InputR(i) + RR * Ef(Nz); % Backward field at z=L with right mirror

    Ef(2:Nz) = fsync * Ef(1:Nz-1); % Forward field propagation
    Er(1:Nz-1) = fsync * Er(2:Nz); % Backward field propagation

    OutputR(i) = Ef(Nz) * (1 - RR); % Transmitted field at z=L
    OutputL(i) = Er(1) * (1 - RL); % Transmitted field at z=0

    % Plot updated fields periodically
    if mod(i, plotN) == 0
        subplot(3, 1, 1)
        plot(z * 10000, real(Ef), 'r'); hold on
        plot(z * 10000, imag(Ef), 'r--'); hold off % Plot imaginary part
        xlim(XL * 1e4)
        ylim(YL)
        xlabel('z(\mum)')
        ylabel('E_f')
        legend('\Re', '\Im')
        hold off

        subplot(3, 1, 2)
        plot(z * 10000, real(Er), 'b'); hold on
        plot(z * 10000, imag(Er), 'b--'); hold off % Plot imaginary part
        xlim(XL * 1e4)
        ylim(YL)
        xlabel('z(\mum)')
        ylabel('E_r')
        legend('\Re', '\Im')
        hold off

        subplot(3, 1, 3)
        plot(time * 1e12, real(InputL), 'r'); hold on
        plot(time * 1e12, real(OutputR), 'g'); % Right output
        plot(time * 1e12, real(InputR), 'b'); % Right input
        plot(time * 1e12, real(OutputL), 'm'); % Left output
        xlim([0, Nt * dt * 1e12])
        ylim(YL)
        xlabel('time(ps)')
        ylabel('E')
        legend('Left Input', 'Right Output', 'Right Input', 'Left Output', 'Location', 'east')
        hold off

        pause(0.01) % Pause for visualization
    end
end
