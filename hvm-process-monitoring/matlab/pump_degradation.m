%% pump_degradation.m
% Vacuum pump wear model and failure forecast for chamber ETCH04.
%
% Model (same as src/reliability.py):
%   dD/dt = k * u(t) * (1 + D/Dc)     damage grows faster as the pump wears
%   V(t)  = V0 + D(t)                 measured vibration velocity (mm/s)
%   u(t)  = 0.85 + 0.15*sin(2*pi*t)   pump duty cycle (t in days)
%
% Steps: fit k, Dc and V0 to the vibration history up to day 40 with
% fminsearch (each evaluation integrates the ODE with ode45), forecast
% forward, and report when vibration reaches the trip level.
%
% Needs base MATLAB only (no toolboxes). Run run_pipeline.py first: it writes
% outputs/pump_ETCH04_first_life.csv, which this script reads.

clear; clc; close all;

here = fileparts(mfilename('fullpath'));
data = readtable(fullfile(here, '..', 'outputs', 'pump_ETCH04_first_life.csv'));
tAll = data.day;
vAll = data.pump_vibration_mm_s;

forecastDay = 40;      % "today": only data up to here is used for the fit
vibTrip     = 7.1;     % protective trip limit (mm/s)
vibAlert    = 4.5;     % alert level (mm/s)
horizon     = 90;      % planning horizon (days)

inTrain = tAll <= forecastDay;
t = tAll(inTrain);
v = vAll(inTrain);

duty = @(tt) 0.85 + 0.15 * sin(2 * pi * tt);

%% Fit with fminsearch. Parameters: p = [log(k), log(Dc), V0] keeps k, Dc > 0
p0   = [log(0.01), log(5), median(v(1:48))];
sse  = @(p) sum((simulateWear(p, t, duty) - v).^2);
opts = optimset('MaxFunEvals', 4000, 'MaxIter', 4000, 'TolX', 1e-8, 'TolFun', 1e-10);
pHat = fminsearch(sse, p0, opts);

k    = exp(pHat(1));
Dc   = exp(pHat(2));
V0   = pHat(3);
rmse = sqrt(sse(pHat) / numel(t));

% The interlock trips on the first single reading over 7.1 mm/s, so the pump is
% at risk once the mean vibration is within two noise SDs of the limit.
tripLevel = vibTrip - 2 * rmse;

%% Forecast
tFull = (0:1/24:forecastDay + horizon)';
vFull = simulateWear(pHat, tFull, duty);

future   = tFull > t(end);
iTrip    = find(future & vFull >= tripLevel, 1);
iAlert   = find(vFull >= vibAlert, 1);
tripDay  = tFull(iTrip);
alertDay = tFull(iAlert);

failRow   = find(vAll >= vibTrip, 1);
actualDay = NaN;
if ~isempty(failRow), actualDay = tAll(failRow); end

fprintf('Fitted parameters:  k = %.4f /day,  Dc = %.3f mm/s,  V0 = %.3f mm/s\n', k, Dc, V0);
fprintf('Fit RMSE: %.3f mm/s   -> trip level used: %.2f mm/s\n', rmse, tripLevel);
fprintf('Alert level (4.5 mm/s) reached on day %.1f\n', alertDay);
fprintf('Forecast made on day %d predicts failure on day %.2f\n', forecastDay, tripDay);
if ~isnan(actualDay)
    fprintf('Actual failure in the data: day %.2f  (error %+.2f days)\n', actualDay, tripDay - actualDay);
end
fprintf('Recommended PM: schedule pump replacement by day %d\n', floor(tripDay - 2));

%% Plot
figure('Color', 'w', 'Position', [100 100 900 420]);
plot(tAll, vAll, 'Color', [0.72 0.83 0.96]); hold on;
plot(t, v, 'Color', [0.38 0.31 0.84]);
plot(tFull(tFull <= tripDay + 2), vFull(tFull <= tripDay + 2), 'k--', 'LineWidth', 1.6);
yline(vibAlert, ':', 'alert 4.5 mm/s', 'Color', [0.79 0.52 0]);
yline(vibTrip,  ':', 'trip 7.1 mm/s',  'Color', [0.82 0.23 0.23]);
xline(forecastDay, '--', 'forecast made', 'Color', [0.54 0.53 0.51]);
xline(tripDay, '-', sprintf('predicted failure day %.1f', tripDay), 'Color', [0 0 0]);
grid on; box off;
xlabel('day'); ylabel('pump vibration (mm/s)');
title('ETCH04 dry pump: ODE wear model (ode45 + fminsearch)');
legend({'measured', 'used for fit', 'model + forecast'}, 'Location', 'northwest');
exportgraphics(gcf, fullfile(here, 'pump_forecast_matlab.png'), 'Resolution', 150);

%% Local function (must sit at the end of the script)
function vib = simulateWear(p, tq, duty)
    k  = exp(p(1));
    Dc = exp(p(2));
    V0 = p(3);
    tq = tq(:);
    addZero = tq(1) > 0;
    if addZero, tq = [0; tq]; end
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    [~, D] = ode45(@(tt, d) k * duty(tt) * (1 + d / Dc), tq, 0, opts);
    if addZero, D = D(2:end); end
    vib = V0 + D;
end
