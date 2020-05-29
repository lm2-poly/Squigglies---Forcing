global Rc SR Uc Ay lambda
SR = 3.5; FR = 3; LR = 5;

fc = 1; Rc = 1; Uc = 2*pi*fc*Rc;
Vp = Uc/SR; fy = fc*FR; Ay = Rc*LR;
                 
lambda = Uc/fy;

tspan = [0 1e3]; Y0 = [1 0 pi/2];
options = odeset('RelTol', 1e-6, 'AbsTol', [1e-5 1e-5 1e-5], 'refine', 5);
[t,Y] = ode23(@forcing, tspan, Y0, options);
Yc = Ay*sin(2*pi*Uc*t/lambda);

X1 = Y(:,1).*cos(Y(:,2));
Y1 = Y(:,1).*sin(Y(:,2)) + Yc;
X2 = X1 + (Uc/SR).*(t(end) - t);
Y2 = Y1;

% Extra
fun = @(t) Vp*sqrt(1+(2*pi*Ay*fy/Vp*cos(2*pi*fy*t)).^2);
lambda2 = integral(fun,0,1/fy);
if LR*FR > pi/2
status = 'WARNING: Uc < Uw'; else       % Coiling speed lower than average wave speed
status = '-'; end
if lambda2 > lambda
status2 = 'WARNING: E < lambda'; else   % Extruded material within a wave lower than the own wave arc length
status2 = '-'; end

%% Plotting
txt = cell(3,1); txt{1} = ['SR = ', num2str(SR)]; txt{2} = ['FR = ', num2str(FR)];
                 txt{3} = ['LR = ', num2str(LR)];
                 
%fname = ['SR', num2str(SR), 'FR', num2str(FR), 'LR', num2str(LR)];
%mkdir('Plots',fname)
%fpath = fullfile(pwd,'Plots',fname);
figure(1)
plot(X1/Rc,Y1/Rc)
title('Sinusoidal movement - Plot of contact point trace')
xlabel('X(s) (mm)', 'FontWeight', 'Bold')
ylabel('Y(s) (mm)', 'FontWeight', 'Bold')
axis([min(X1) max(X1), min(Y1) max(Y1)])
text(min(X1)+0.04*(max(X1)-min(X1)), max(Y1)-0.10*(max(Y1)-min(Y1)), txt, 'FontWeight', 'Bold', 'FontSize', 9)
grid
%saveas(figure(1),fullfile(fpath,'plot1.fig'))

figure(2)
plot(X2/Rc,Y2/Rc,'-r')
title('Sinusoidal movement - Plot of deposited trace')
xlabel('X(s) (mm)', 'FontWeight', 'Bold')
ylabel('Y(s) (mm)', 'FontWeight', 'Bold')
axis([min(X2) max(X2), min(Y2) max(Y2)])
text(min(X2)+0.04*(max(X2)-min(X2)), max(Y2)-0.10*(max(Y2)-min(Y2)), txt, 'FontWeight', 'Bold', 'FontSize', 9)
grid
%saveas(figure(2),fullfile(fpath,'plot2.fig'))

%%
function dY = forcing(t,Y)
global Rc SR Uc Ay lambda
dY = zeros(3,1);
b = 0.715;
A = b^2 * cos(Y(3)-Y(2)) / (1 - b*cos(Y(3)-Y(2)));
kappa = sqrt(Y(1)/Rc^3) * (1 + A*Y(1)/Rc) * sin(Y(3)-Y(2));
dY(1) = Uc*(cos(Y(3)-Y(2)) + (1/SR)*cos(Y(2)) - 2*pi*Ay/lambda*cos(2*pi*Uc*t/lambda)*sin(Y(2)));
dY(2) = Uc*(sin(Y(3)-Y(2)) - (1/SR)*sin(Y(2)) - 2*pi*Ay/lambda*cos(2*pi*Uc*t/lambda)*cos(Y(2))) / Y(1);
dY(3) = Uc*kappa;
end