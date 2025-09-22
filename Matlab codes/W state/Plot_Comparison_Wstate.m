clear all;
close all;
clc;
syms al

kt_vals = 0:0.1:1;
q1=0.1;
q2=0.9;
for jj = 1:length(kt_vals)
        kt = kt_vals(jj);
        g_WM(jj)=exp(-3*kt)*(q1 - 1)^3*(q1 + exp(-kt) - q1*exp(-kt) - 2);
        F_WM(jj)=   -(exp(-4*kt)*(q1 - 1)^3*(4*q1 + 19*exp(kt) - 4*q1*exp(kt) - 4))/15;

        g_WM2(jj)=exp(-3*kt)*(q2 - 1)^3*(q2 + exp(-kt) - q2*exp(-kt) - 2);
        F_WM2(jj)=   -(exp(-4*kt)*(q2 - 1)^3*(4*q2 + 19*exp(kt) - 4*q2*exp(kt) - 4))/15;

        F_NP(jj)=(4*abs(exp(-kt) - 1))/15 + exp(-kt);
        g_s(jj)=exp(-3*kt);
        f_s(jj)=1;
        g_NP(jj)=1;
end

figure(1);
L1 = plot(kt_vals, F_WM, 'm--', 'LineWidth', 1.5); hold on; % Magenta color
L2 =plot(kt_vals, F_WM2, 'Color', [0 0.5 0], 'LineStyle', '--', 'LineWidth', 2); hold on
L3 = plot(kt_vals, f_s, 'b-', 'LineWidth', 1.5);  % Blue color
L4 = plot(kt_vals, F_NP, 'r:', 'LineWidth', 2);  % Black color
axis tight
grid on

legend([L1, L2 ,L3, L4], {'WMR_{0.1}', 'WMR_{0.9}', 'ES=EAM', 'Unprotected'}, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

xlim([0 1]);
xticks(0:0.2:1);
xlabel('$\lambda t$','Interpreter', 'latex',  'FontWeight', 'bold', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Average teleportation fidelity', 'FontWeight', 'bold', 'FontName', 'Times New Roman', 'FontSize', 12);
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';


figure(2);
L5=plot(kt_vals,g_WM, 'm--', 'LineWidth', 1.5); hold on
L6=plot(kt_vals, g_WM2, 'Color', [0 0.5 0], 'LineStyle', '--', 'LineWidth', 2); hold on
L7=plot(kt_vals,g_s, 'b-', 'LineWidth', 1.5); hold on
L8=plot(kt_vals,g_NP, 'r:', 'LineWidth', 2);
axis tight
grid on

legend([L5, L6, L7, L8], {'WMR_{0.1}', 'WMR_{0.9}', 'ES=EAM','Unprotected'}, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

xlim([0 1]);
xticks(0:0.2:1);
xlabel('$\lambda t$','Interpreter', 'latex',  'FontWeight', 'bold', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Total teleportation success probability', 'FontWeight', 'bold', 'FontName', 'Times New Roman', 'FontSize', 12);
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
