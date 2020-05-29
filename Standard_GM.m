global Rc SR Uc
SR = 3.4;

fc = 1; Rc = 1; Uc = 2*pi*fc*Rc;
Vp = Uc/SR;

tspan = [0 1e3]; Y0 = [1 0 pi/2];
options = odeset('RelTol', 1e-6, 'AbsTol', [1e-5 1e-5 1e-5], 'refine', 5);
[t,Y] = ode23(@noforcing, tspan, Y0, options);

X1 = Y(:,1).*cos(Y(:,2));
Y1 = Y(:,1).*sin(Y(:,2));
X2 = X1 + (Uc/SR).*(t(end) - t);
Y2 = Y1;


%% Plotting
txt = cell(1,1); txt{1} = ['SR = ', num2str(SR)];
n = round(0.7*length(Y)):length(Y); m = round(0.95*length(Y)):length(Y);

figure(1)
plot(X1(n)/Rc,Y1(n)/Rc)
title('Sinusoidal movement - Plot of contact point trace')
xlabel('X(s)', 'FontWeight', 'Bold')
ylabel('Y(s)', 'FontWeight', 'Bold')

figure(2)
plot(X2(m)/Rc,Y2(m)/Rc)
title('Sinusoidal movement - Plot of deposited trace')
xlabel('X(s)', 'FontWeight', 'Bold')
ylabel('Y(s)', 'FontWeight', 'Bold')

%%
function dY = noforcing(t,Y)
global Rc SR Uc
dY = zeros(3,1);
b = 0.715;
A = b^2 * cos(Y(3)-Y(2)) / (1 - b*cos(Y(3)-Y(2)));
kappa = sqrt(Y(1)/Rc^3) * (1 + A*Y(1)/Rc) * sin(Y(3)-Y(2));
dY(1) = Uc*(cos(Y(3)-Y(2)) + (1/SR)*cos(Y(2)));
dY(2) = Uc*(sin(Y(3)-Y(2)) - (1/SR)*sin(Y(2))) / Y(1);
dY(3) = Uc*kappa;
end