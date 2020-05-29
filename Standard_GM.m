global Rc SR Uc
SR = 3.4;                                                                         % Speed ratio

fc = 1; Rc = 1; Uc = 2*pi*fc*Rc;                                                  % Steady coiling frequency, radius and extruding speed
Vp = Uc/SR;                                                                       % Belt speed

tspan = [0 1e3]; Y0 = [1 0 pi/2];                                                 % Temporal integration interval and initial conditions
options = odeset('RelTol', 1e-6, 'AbsTol', [1e-5 1e-5 1e-5], 'refine', 5);        % ODE options
[t,Y] = ode23(@noforcing, tspan, Y0, options);                                    % ODE numerical resolution

X1 = Y(:,1).*cos(Y(:,2)); X2 = X1 + (Uc/SR).*(t(end) - t);                        % X-coordinate contact point trace / deposited trace                      
Y1 = Y(:,1).*sin(Y(:,2)); Y2 = Y1;                                                % Y-coordinate contact point trace / deposited trace


%% ---------------------------------------------------- PLOTTING -----------------------------------------------------------------------------------------------------
txt = cell(1,1); txt{1} = ['SR = ', num2str(SR)];                                 % textbox referencing the speed ratio used
n = round(0.7*length(Y)):length(Y); m = round(0.95*length(Y)):length(Y);          % cutting solution array allows to work on the steady spectrum

figure(1)                                                                         % x- and y- coordinates are normalized over Rc
plot(X1(n)/Rc,Y1(n)/Rc)                                                           % Plot of contact point trace
title('Sinusoidal movement - Plot of contact point trace')
xlabel('y(s)', 'FontWeight', 'Bold'); ylabel('y(s)', 'FontWeight', 'Bold')

figure(2)                                                                         % X- and Y- coordinates are normalized over Rc
plot(X2(m)/Rc,Y2(m)/Rc)                                                           % Plot of deposited trace
title('Sinusoidal movement - Plot of deposited trace')
xlabel('X(s)', 'FontWeight', 'Bold'); ylabel('Y(s)', 'FontWeight', 'Bold')

%% ---------------------------------------------------- CALLED FUNCTIONS -----------------------------------------------------------------------------------------------------
function dY = noforcing(t,Y)                                                      % Y(1) = r, Y(2) = phi, Y(3) = theta                                 
global Rc SR Uc
dY = zeros(3,1);
b = 0.715;
A = b^2 * cos(Y(3)-Y(2)) / (1 - b*cos(Y(3)-Y(2)));
kappa = sqrt(Y(1)/Rc^3) * (1 + A*Y(1)/Rc) * sin(Y(3)-Y(2));
dY(1) = Uc*(cos(Y(3)-Y(2)) + (1/SR)*cos(Y(2)));
dY(2) = Uc*(sin(Y(3)-Y(2)) - (1/SR)*sin(Y(2))) / Y(1);
dY(3) = Uc*kappa;
end
