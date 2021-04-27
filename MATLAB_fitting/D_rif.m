function [fitresult, gof] = D_rif(time, meanC_1, meanC_2, meanC_3, meanC_4, meanC_5, meanC_6, meanC_7)
%CREATEFITS1(TIME,MEANC_1,MEANC_2,MEANC_3,MEANC_4,MEANC_5,MEANC_6,MEANC_7)
%  Create fits.
%
%  Data for 'C=0' fit:
%      X Input : time
%      Y Output: meanC_1
%  Data for 'C=0.0625' fit:
%      X Input : time
%      Y Output: meanC_2
%  Data for 'C=0.125' fit:
%      X Input : time
%      Y Output: meanC_3
%  Data for 'C=0.250' fit:
%      X Input : time
%      Y Output: meanC_4
%  Data for 'C=0.5' fit:
%      X Input : time
%      Y Output: meanC_5
%  Data for 'C=1' fit:
%      X Input : time
%      Y Output: meanC_6
%  Data for 'C=2' fit:
%      X Input : time
%      Y Output: meanC_7
%  Output:
%      fitresult : a cell-array of fit objects representing the fits.
%      gof : structure array with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 04-Oct-2018 03:56:49

%% Initialization.

% Initialize arrays to store fits and goodness-of-fit.
fitresult = cell( 7, 1 );
gof = struct( 'sse', cell( 7, 1 ), ...
    'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );

%% Fit: 'C=0'.
[xData, yData] = prepareCurveData( time, meanC_1 );

% Set up fittype and options.
ft = fittype( '(k1./(1 + (k1-n0)./n0*exp(-r1.*(x-xc))))*heaviside(xc-x) + (k2./(1 + (k2-n0)./n0*exp(-r2.*(x-xc))))*heaviside(x-xc)', 'independent', 'x', 'dependent', 'y' );
excludedPoints = excludedata( xData, yData, 'Indices', [36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62] );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.8 0.8 0.4 0.111202755293787 0.780252068321138 7.5];
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult{1}, gof(1)] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'C=0' );
plot( fitresult{1}, xData, yData, excludedPoints );
% Label axes
xlabel time
ylabel meanC_1
grid on

%% Fit: 'C=0.0625'.
[xData, yData] = prepareCurveData( time, meanC_2 );

% Set up fittype and options.
ft = fittype( '(k1./(1 + (k1-n0)./n0*exp(-r1.*(x-xc))))*heaviside(xc-x) + (k2./(1 + (k2-n0)./n0*exp(-r2.*(x-xc))))*heaviside(x-xc)', 'independent', 'x', 'dependent', 'y' );
excludedPoints = excludedata( xData, yData, 'Indices', [36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62] );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.4 1 0.4 0.0964545251683886 0.131973292606335 7];
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult{2}, gof(2)] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'C=0.0625' );
plot( fitresult{2}, xData, yData, excludedPoints );
% Label axes
xlabel time
ylabel meanC_2
grid on

%% Fit: 'C=0.125'.
[xData, yData] = prepareCurveData( time, meanC_3 );

% Set up fittype and options.
ft = fittype( '(k1./(1 + (k1-n0)./n0*exp(-r1.*(x-xc))))*heaviside(xc-x) + (k2./(1 + (k2-n0)./n0*exp(-r2.*(x-xc))))*heaviside(x-xc)', 'independent', 'x', 'dependent', 'y' );
excludedPoints = excludedata( xData, yData, 'Indices', [36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62] );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.4 0.575208595078466 0.438744359656398 0.0597795429471558 0.234779913372406 10];
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult{3}, gof(3)] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'C=0.125' );
plot( fitresult{3}, xData, yData, excludedPoints );
% Label axes
xlabel time
ylabel meanC_3
grid on

%% Fit: 'C=0.250'.
[xData, yData] = prepareCurveData( time, meanC_4 );

% Set up fittype and options.
ft = fittype( '(k1./(1 + (k1-n0)./n0*exp(-r1.*(x-xc))))*heaviside(xc-x) + (k2./(1 + (k2-n0)./n0*exp(-r2.*(x-xc))))*heaviside(x-xc)', 'independent', 'x', 'dependent', 'y' );
excludedPoints = excludedata( xData, yData, 'Indices', [1 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62] );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0 0 0 0];
opts.StartPoint = [0.5 1 0.4456 1 0.168990029462704 7];
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult{4}, gof(4)] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'C=0.250' );
plot( fitresult{4}, xData, yData, excludedPoints );
% Label axes
xlabel time
ylabel meanC_4
grid on

%% Fit: 'C=0.5'.
[xData, yData] = prepareCurveData( time, meanC_5 );

% Set up fittype and options.
ft = fittype( '(k1./(1 + (k1-n0)./n0*exp(-r1.*(x-xc))))*heaviside(xc-x) + (k2./(1 + (k2-n0)./n0*exp(-r2.*(x-xc))))*heaviside(x-xc)', 'independent', 'x', 'dependent', 'y' );
excludedPoints = excludedata( xData, yData, 'Indices', [36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62] );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0 0 0 0];
opts.StartPoint = [0.73172238565867 0.647745963136307 0 0.450923706430945 0.547008892286345 10];
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult{5}, gof(5)] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'C=0.5' );
plot( fitresult{5}, xData, yData, excludedPoints );
% Label axes
xlabel time
ylabel meanC_5
grid on

%% Fit: 'C=1'.
[xData, yData] = prepareCurveData( time, meanC_6 );

% Set up fittype and options.
ft = fittype( '(k1./(1 + (k1-n0)./n0*exp(-r1.*(x-xc))))*heaviside(xc-x) + (k2./(1 + (k2-n0)./n0*exp(-r2.*(x-xc))))*heaviside(x-xc)', 'independent', 'x', 'dependent', 'y' );
excludedPoints = excludedata( xData, yData, 'Indices', [36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62] );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.5 0.188955015032545 0.5853 0.686775433365315 0.18351115573727 10];
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult{6}, gof(6)] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'C=1' );
plot( fitresult{6}, xData, yData, excludedPoints );
% Label axes
xlabel time
ylabel meanC_6
grid on

%% Fit: 'C=2'.
[xData, yData] = prepareCurveData( time, meanC_7 );

% Set up fittype and options.
ft = fittype( '(k1./(1 + (k1-n0)./n0*exp(-r1.*(x-xc))))*heaviside(xc-x) + (k2./(1 + (k2-n0)./n0*exp(-r2.*(x-xc))))*heaviside(x-xc)', 'independent', 'x', 'dependent', 'y' );
excludedPoints = excludedata( xData, yData, 'Indices', [36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62] );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.62561856072969 0.780227435151377 0.890903252535798 0.0811257688657853 0.92938597096873 10];
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult{7}, gof(7)] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'C=2' );
plot( fitresult{7}, xData, yData, excludedPoints );
% Label axes
xlabel time
ylabel meanC_7
grid on


