global Rc SR Uc Ay lambda
SR = 3.5; FR = 3; LR = 5;                                                            % Speed, frequency and length ratios

fc = 1; Rc = 1; Uc = 2*pi*fc*Rc;                                                     % Steady coiling frequency, radius and extruding speed
Vp = Uc/SR; fy = fc*FR; Ay = Rc*LR;                                                  % Belt speed, forcing frequency and forcing amplitude
                                                                                     % Note: FR is with respect to fc (See "Complete.m" for FR with respect to pattern frequency)
lambda = Uc/fy;                                                                      % lambda is the total amount of material extruded within a perturbation cycle

tspan = [0 1e3]; Y0 = [1 0 pi/2];                                                    % Temporal integration interval and initial conditions
options = odeset('RelTol', 1e-6, 'AbsTol', [1e-5 1e-5 1e-5], 'refine', 5);           % ODE options
[t,Y] = ode23(@forcing, tspan, Y0, options);                                         % ODE numerical resolution - FORCING
Yc = Ay*sin(2*pi*Uc*t/lambda);                                                       % Sinusoidal forcing function

X1 = Y(:,1).*cos(Y(:,2)); X2 = X1 + (Uc/SR).*(t(end) - t);                           % X-coordinate contact point trace / deposited trace
Y1 = Y(:,1).*sin(Y(:,2)) + Yc; Y2 = Y1;                                              % Y-coordinate contact point trace / deposited trace

%% ---------------------------------------------------- OTHER CONSIDERATIONS -----------------------------------------------------------------------------------------------------
fun = @(t) Vp*sqrt(1+(2*pi*Ay*fy/Vp*cos(2*pi*fy*t)).^2);
lambda2 = integral(fun,0,1/fy);                                                      % Regular sine wave arc length
if LR*FR > pi/2
status = 'WARNING: Uc < Uw'; else                                                    % Coiling speed is lower than average wave speed
status = '-'; end
if lambda2 > lambda
status2 = 'WARNING: E < lambda'; else                                                % Extruded material within a wave is lower than the own wave arc length
status2 = '-'; end

%% ---------------------------------------------------- PLOTTING -----------------------------------------------------------------------------------------------------
n = round(length(Y0)*0.5):length(Y0); m = round(length(Y0)*0.9):length(Y0);          % cutting solution array allows to work on the steady spectrum, if any.

figure(1)
plot(X1(n)/Rc,Y1(n)/Rc)                                                              % x- and y-coordinates are normalized over Rc
title('Sinusoidal movement - Plot of contact point trace')
xlabel('x(s)', 'FontWeight', 'Bold')
ylabel('y(s)', 'FontWeight', 'Bold')

figure(2)
plot(X2(m)/Rc(m),Y2/Rc,'-r')                                                         % X- and Y-coordinates are normalized over Rc
title('Sinusoidal movement - Plot of deposited trace')
xlabel('X(s)', 'FontWeight', 'Bold')
ylabel('Y(s)', 'FontWeight', 'Bold')

%% ---------------------------------------------------- CALLED FUNCTIONS -----------------------------------------------------------------------------------------------------
function dY = forcing(t,Y)                                                           % Y(1) = r, Y(2) = phi, Y(3) = theta
global Rc SR Uc Ay lambda
dY = zeros(3,1);
b = 0.715;
A = b^2 * cos(Y(3)-Y(2)) / (1 - b*cos(Y(3)-Y(2)));
kappa = sqrt(Y(1)/Rc^3) * (1 + A*Y(1)/Rc) * sin(Y(3)-Y(2));
dY(1) = Uc*(cos(Y(3)-Y(2)) + (1/SR)*cos(Y(2)) - 2*pi*Ay/lambda*cos(2*pi*Uc*t/lambda)*sin(Y(2)));
dY(2) = Uc*(sin(Y(3)-Y(2)) - (1/SR)*sin(Y(2)) - 2*pi*Ay/lambda*cos(2*pi*Uc*t/lambda)*cos(Y(2))) / Y(1);
dY(3) = Uc*kappa;
end
