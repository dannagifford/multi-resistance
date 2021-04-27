%% read mat file

clear
close all

load('blanked-growth-curve-data.mat')

%%

% specify antibiotic and strain
% indexA specifies the antibiotic from antibiotic_array
% indexS specifies the strain from strain_array
% e.g. for strain A under rifampicin set indexA = 1, indexS = 2;

antibiotic_array = {'rifampicin', 'nalidixic acid', 'combination'}; indexA = 1;
strain_array = {'S', 'A', 'B', 'D'}; indexS = 2;

strain = strain_array{indexS};
antibiotic = antibiotic_array{indexA};

C_array = [0 0.0625 0.125 0.25 0.5 1 2];

% contruct matrix with the mean values of every concentration C

single_time = 62;

meanC = ones(length(C_array),single_time);

figure(1)
hold on
for i=1:1:length(C_array)
    [time, meanC(i, :)] = extractMean(blankedgrowthcurvedata, antibiotic, C_array(i), strain);
    plot(time, meanC(i, :), 'LineWidth', 2)
end

box on
grid on
set(gca, 'Fontsize', 30, 'LineWidth', 3, 'TickLabelInterpreter', 'latex')
xlabel('time (hours)', 'Interpreter', 'latex')
ylabel('OD', 'Interpreter', 'latex')
l=legend({'C=0', 'C=0.625','C=1.25', 'C=2.5', 'C=5', 'C=10', 'C=20'}, 'Interpreter', 'latex');
l.Location = 'best';
axis tight

meanC_1 = transpose(meanC(1, :));
meanC_2 = transpose(meanC(2, :));
meanC_3 = transpose(meanC(3, :));
meanC_4 = transpose(meanC(4, :));
meanC_5 = transpose(meanC(5, :));
meanC_6 = transpose(meanC(6, :));
meanC_7 = transpose(meanC(7, :));

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gcf,'rend','painters','pos',[50 50 900 600])

%%
% change here the strain and treatment considered 
% e.g. A_rif(...), S_nal(...), D_both(...), etc

[fitresult, gof] = A_rif(time, meanC_1, meanC_2, meanC_3, meanC_4, meanC_5, meanC_6, meanC_7);

%% R^2 of fitting

R2=extractfield(gof,'rsquare');

%% generate k, r, and n0 arrays with the values from the fitting

k1 = zeros(7,1); k2 = zeros(7,1);
n0 = zeros(7,1); xc = zeros(7,1);
r1 = zeros(7,1); r2 = zeros(7,1);
k1_error_U = zeros(7,1); k2_error_U = zeros(7,1);
k1_error_L = zeros(7,1); k2_error_L = zeros(7,1);
r1_error_U = zeros(7,1); r2_error_U = zeros(7,1);
r1_error_L = zeros(7,1); r2_error_L = zeros(7,1);
n0_error_L = zeros(7,1); n0_error_U = zeros(7,1);

for i=1:1:7
    fitting_coeff = coeffvalues(fitresult{i});
    ci = confint(fitresult{i}); %confidence intervals
    k1(i) = fitting_coeff(1); k2(i) = fitting_coeff(2);
    n0(i) = fitting_coeff(3); xc(i) = fitting_coeff(6);
    r1(i) = fitting_coeff(4); r2(i) = fitting_coeff(5);
    k1_error_L(i) = ci(1,1); k1_error_U(i) = ci(2,1);
    k2_error_L(i) = ci(1,2); k2_error_U(i) = ci(2,2);
    n0_error_L(i) = ci(1,3); n0_error_U(i) = ci(2,3);
    r1_error_L(i) = ci(1,4); r1_error_U(i) = ci(2,4);
    r2_error_L(i) = ci(1,5); r2_error_U(i) = ci(2,5);
end


%% plot k 

figure
hold on
errorbar(1:1:7,k1,k1-k1_error_L, k1_error_U-k1,'-*','LineWidth', 2)
set(gca,'Fontsize',17,'LineWidth',1.5)
xlabel('day')
ylabel('k_1')
title(strcat('strain:',{' '},strain, ', antibiotic:',{' '}, antibiotic))
grid on

figure
hold on
errorbar(1:1:7,k2,k2-k2_error_L, k2_error_U-k2,'-*','LineWidth', 2)
set(gca,'Fontsize',17,'LineWidth',1.5)
xlabel('day')
ylabel('k_2')
title(strcat('strain:',{' '},strain, ', antibiotic:',{' '}, antibiotic))
grid on


%% plot r

figure
hold on
errorbar(1:1:7,r1,r1-r1_error_L, r1_error_U-r1,'-*','LineWidth', 2)
set(gca,'Fontsize',17,'LineWidth',1.5)
xlabel('day')
ylabel('r_1')
title(strcat('strain:',{' '},strain, ', antibiotic:',{' '}, antibiotic))
grid on

figure
hold on
errorbar(1:1:7,r2,r2-r2_error_L, r2_error_U-r2,'-*','LineWidth', 2)
set(gca,'Fontsize',17,'LineWidth',1.5)
xlabel('day')
ylabel('r_2')
title(strcat('strain:',{' '},strain, ', antibiotic:',{' '}, antibiotic))
grid on


%% plot n0

figure
hold on
errorbar(1:1:7,n0,n0-n0_error_L, n0_error_U-n0,'-*','LineWidth', 2)
set(gca,'Fontsize',17,'LineWidth',1.5)
xlabel('day')
ylabel('n_0')
title(strcat('strain:',{' '},strain, ', antibiotic:',{' '}, antibiotic))
grid on

