clear all;
close all;
clc;
syms al

kt_vals = 0:0.1:1;
q=0.5;
for jj = 1:length(kt_vals)
        kt = kt_vals(jj);
        g_WM(jj) =(exp(-4*kt)*(q - 1)^2*(5*exp(2*kt) - 2*q - 4*exp(kt) - 4*q*exp(2*kt) - 2*q^2*exp(kt) + q^2*exp(2*kt) + 6*q*exp(kt) + q^2 + 1))/2;
        F_WM(jj)=(49*exp(2*kt) - 22*q - 30*exp(kt) - 30*q*exp(2*kt) - 22*q^2*exp(kt) + 11*q^2*exp(2*kt) + 52*q*exp(kt) + 11*q^2 + 11)/(15*(5*exp(2*kt) - 2*q - 4*exp(kt) - 4*q*exp(2*kt) - 2*q^2*exp(kt) + q^2*exp(2*kt) + 6*q*exp(kt) + q^2 + 1));        
        F_NP(jj)=(7*exp(-2*kt))/15 - exp(-kt)/5 + 11/15;
        g_s(jj)=exp(-2*kt);
        f_s(jj)=1;
        g_NP(jj)=1;
end

figure(1);
L1 = plot(kt_vals, F_WM, 'm--', 'LineWidth', 1.5); hold on; % Magenta color
L2 = plot(kt_vals, f_s, 'b-', 'LineWidth', 1.5);  % Blue color
L3 = plot(kt_vals, F_NP, 'r:', 'LineWidth', 2);  % Black color
axis tight
grid on

legend([L1, L2 ,L3], {'WMR', 'ES=EAM', 'Unprotected'}, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

xlim([0 1]);
xticks(0:0.2:1);
xlabel('$\lambda t$','Interpreter', 'latex',  'FontWeight', 'bold', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Average teleportation fidelity', 'FontWeight', 'bold', 'FontName', 'Times New Roman', 'FontSize', 12);
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';


figure(2);
L3=plot(kt_vals,g_WM, 'm--', 'LineWidth', 1.5); hold on
L4=plot(kt_vals,g_s, 'b-', 'LineWidth', 1.5); hold on
L5=plot(kt_vals,g_NP, 'r:', 'LineWidth', 2);
axis tight
grid on

legend([L3, L4, L5], {'WMR', 'ES=EAM','Unprotected'}, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

xlim([0 1]);
xticks(0:0.2:1);
xlabel('$\lambda t$','Interpreter', 'latex',  'FontWeight', 'bold', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Total teleportation success probability', 'FontWeight', 'bold', 'FontName', 'Times New Roman', 'FontSize', 12);
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
