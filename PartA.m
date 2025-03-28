clc;  % Clears the command window
clear;  % Clears all variables from the workspace
%close all;  % Closes all open figures

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Energy Efficient (EE)Motor parameters')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
 
f_EE=50;               %Supply frequency [Hz]
p_EE=4;                %Number of poles
V1_EE=380/sqrt(3);     %Supply voltage [phase]
R1_EE=1.5;             %Stator winding resistance [ohms/phase]
X1_EE=3.642;           %Stator winding leakage reactance [ohms/phase]
Xm_EE=72.252;          %Stator winding magnetising reactance [ohms/phase]
X2p_EE=3.642;          %Rotor winding leakage reactance reffered to stator [ohms/phase]
R2p_EE=1.994;          %Rotor winding resistance reffered to stator [ohms/phase]

fprintf('EE Motor\n');
fprintf('\n');

fprintf('f=%f\n',f_EE);
fprintf('p=%f\n',p_EE);
fprintf('V1=%f\n',V1_EE);
fprintf('R1=%f\n',R1_EE);
fprintf('X1=%f\n',X1_EE);
fprintf('Xm=%f\n',Xm_EE);
fprintf('X2p=%f\n',X2p_EE);
fprintf('R2p=%f\n',R2p_EE);

f_SE=50;               %Supply frequency [Hz]
p_SE=4;                %Number of poles
V1_SE=380/sqrt(3);     %Supply voltage [phase]
R1_SE=2.087;             %Stator winding resistance [ohms/phase]
X1_SE=4.274;           %Stator winding leakage reactance [ohms/phase]
Xm_SE=66.56;          %Stator winding magnetising reactance [ohms/phase]
X2p_SE=4.2742;          %Rotor winding leakage reactance reffered to stator [ohms/phase]
R2p_SE=2.122;          %Rotor winding resistance reffered to stator [ohms/phase]

fprintf('\n');
fprintf('SE Motor\n');
fprintf('\n');

fprintf('f=%f\n',f_SE);
fprintf('p=%f\n',p_SE);
fprintf('V1=%f\n',V1_SE);
fprintf('R1=%f\n',R1_SE);
fprintf('X1=%f\n',X1_SE);
fprintf('Xm=%f\n',Xm_SE);
fprintf('X2p=%f\n',X2p_SE);
fprintf('R2p=%f\n',R2p_SE);

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Question 1:Thevenin Equiv Cct Parameters for EE and SE Motor:')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

Vth_EE=Xm_EE/sqrt(R1_EE^2+(X1_EE+Xm_EE)^2)*V1_EE;         %Thevenin equiv voltage source [V] (Equ 5.45 - Sen)
Zth_EE=1i*Xm_EE*(R1_EE+1i*X1_EE)/(R1_EE+1i*(X1_EE+Xm_EE));   %Thevenin equiv impedance
Rth_EE=real(Zth_EE);                          %Thevenin equiv resistance [ohms]
Xth_EE=imag(Zth_EE);                          %Thevenin equiv reactance [ohms]

fprintf('EE Motor\n');
fprintf('\n');

fprintf('Vth=%f\n',Vth_EE);
fprintf('Rth=%f\n',Rth_EE);
fprintf('Xth=%f\n',Xth_EE);

Vth_SE=Xm_SE/sqrt(R1_SE^2+(X1_SE+Xm_SE)^2)*V1_SE;         %Thevenin equiv voltage source [V] (Equ 5.45 - Sen)
Zth_SE=1i*Xm_SE*(R1_SE+1i*X1_SE)/(R1_SE+1i*(X1_SE+Xm_SE));   %Thevenin equiv impedance
Rth_SE=real(Zth_SE);                          %Thevenin equiv resistance [ohms]
Xth_SE=imag(Zth_SE);                          %Thevenin equiv reactance [ohms]

fprintf('\n');
fprintf('SE Motor\n');
fprintf('\n');

fprintf('Vth=%f\n',Vth_SE);
fprintf('Rth=%f\n',Rth_SE);
fprintf('Xth=%f\n',Xth_SE);

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('QUESTION 2:Torque versus speed characteristics for EE and SE Motor:')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

ns=120*f_EE/p_EE;         %Synchronous speed [rpm]
ws=2*pi*ns/60;      %Synchronous speed [rad/sec]
s=0.0005:0.0005:1;  %Slip [pu]
n=(1-s)*ns;         %Rotor speed [rpm]
w=2*pi*n/60;        %Rotor speed [rad/sec]
 
Tmech_EE=3/ws*Vth_EE^2./((Rth_EE+R2p_EE./s).^2+(Xth_EE+X2p_EE)^2).*R2p_EE./s;    %Total Tmech = {3*(Equ5.54 - Sen)}
Tmech_SE=3/ws*Vth_SE^2./((Rth_SE+R2p_SE./s).^2+(Xth_SE+X2p_SE)^2).*R2p_SE./s;

% Find starting and maximum torque
% Starting torque: s = 1
T_start_EE = 3/ws*Vth_EE^2./((Rth_EE+R2p_EE./1).^2+(Xth_EE+X2p_EE)^2).*R2p_EE./1;
T_start_SE = 3/ws*Vth_SE^2./((Rth_SE+R2p_SE./1).^2+(Xth_SE+X2p_SE)^2).*R2p_SE./1;

% Max torque
T_max_EE = (3/(2*ws)).*Vth_EE^2./((Rth_EE+(Rth_EE^2 +(Xth_EE+X2p_EE)^2)^0.5));
T_max_SE = (3/(2*ws)).*Vth_SE^2./((Rth_SE+(Rth_SE^2 +(Xth_SE+X2p_SE)^2)^0.5));

st_max_EE = R2p_EE/(Rth_EE^2+(Xth_EE+X2p_EE)^2)^0.5;
n_max_EE = (1 - st_max_EE) * ns;

st_max_SE = R2p_SE/(Rth_SE^2+(Xth_SE+X2p_SE)^2)^0.5;
n_max_SE = (1 - st_max_SE) * ns;

% Plot Torque vs Speed for EE and SE Motors
figure;
subplot(2,2,1),
plot(n, Tmech_EE, 'r', 'LineWidth', 2), hold on;
plot(n, Tmech_SE, 'b', 'LineWidth', 2);
xlabel('Rotor Speed (RPM)'), ylabel('Torque (Nm)'),...
title('Torque vs Speed for EE and SE Motors'), grid on;
legend('EE Motor', 'SE Motor');

fprintf('\n2.a)\n Starting torque for EE Motor: %.4f Nm\n', T_start_EE);
fprintf(' Starting torque for SE Motor: %.4f Nm\n\n', T_start_SE);

disp(' The starting torque will vary with a change in the rotor resistance R2p.');
disp(' The starting torque is also proportional to the square of supply voltage Vth.');

fprintf('\n2.b)\n Maximum torque for EE Motor: %.4f Nm\n', T_max_EE);
fprintf(' Maximum torque for SE Motor: %.4f Nm\n\n', T_max_SE);

disp('The maximum torque will vary with a change in the supply voltage Vth, as it is independent of R2p.');

fprintf('\n2.c) Speed at which maximum torque occurs for EE Motor: %.4f rpm\n', n_max_EE);
fprintf('Speed at which maximum torque occurs for SE Motor: %.4f rpm\n\n', n_max_SE);

disp('The value of the rotor circuit resistance R2p determines the speed at which maximum torque will occur.');


% Extend the code to calculate the stator current vs. speed
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('QUESTION 3: Stator Current vs. Speed Characteristics for EE and SE Motor:')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

% Stator current calculation
Z1_EE = R1_EE+(1i*X1_EE)+1i*Xm_EE*((R2p_EE./s)+1i*X2p_EE)./((R2p_EE./s)+1i*(Xm_EE+X2p_EE)); %(Equ5.65a - Sen)
I1_EE = V1_EE./Z1_EE;  %(Equ5.65c - Sen) 
I1_mag_EE = abs(I1_EE);
I1_phase_EE = angle(I1_EE);

Z1_SE = R1_SE+(1i*X1_SE)+1i*Xm_SE*((R2p_SE./s)+1i*X2p_SE)./((R2p_SE./s)+1i*(Xm_SE+X2p_SE)); %(Equ5.65a - Sen)
I1_SE = V1_SE./Z1_SE;  %(Equ5.65c - Sen) 
I1_mag_SE = abs(I1_SE);
I1_phase_SE = angle(I1_SE);

% Plot Stator Current vs Speed for EE and SE Motors
figure;
subplot(2,2,2),
plot(n, I1_mag_EE, 'r', 'LineWidth', 2), hold on;
plot(n, I1_mag_SE, 'b', 'LineWidth', 2);
xlabel('Rotor Speed [RPM]'), ylabel('Stator Current [A]'),...
title('Stator Current vs Speed for EE and SE Motors'), grid on;
legend('EE Motor', 'SE Motor');

% a) Stator current at start-up (s=1)

Z1_start_EE = R1_EE+(1i*X1_EE)+1i*Xm_EE*((R2p_EE./1)+1i*X2p_EE)./((R2p_EE./1)+1i*(Xm_EE+X2p_EE)); %(Equ5.65a - Sen)
I_start_EE = V1_EE./Z1_start_EE;

Z1_start_SE = R1_SE+(1i*X1_SE)+1i*Xm_SE*((R2p_SE./1)+1i*X2p_SE)./((R2p_SE./1)+1i*(Xm_SE+X2p_SE)); %(Equ5.65a - Sen)
I_start_SE = V1_SE./Z1_start_SE;

fprintf('\n3.a) Stator current at start-up\n');
fprintf('Stator current for EE Motor at start-up: %.4f < %.2f° A\n', abs(I_start_EE), rad2deg(angle(I_start_EE)));
fprintf('Stator current for SE Motor at start-up: %.4f < %.2f° A\n', abs(I_start_SE), rad2deg(angle(I_start_SE)));

fprintf('\n');
disp('At start-up, the stator current is high because the slip is 1, leading to higher rotor resistance and higher current draw.');

% b) Explain stator current change under no-load and full-load conditions
disp('b) The stator current at start-up will be highest due to the high slip (s=1). Under no-load conditions, the slip is small (close to 0), and the current decreases significantly, as most of the current is used to magnetize the machine. Under full-load conditions, the slip increases, and the stator current increases again due to higher losses in the rotor resistance.');

% c) Stator current at maximum torque
Z1_tmax_EE = R1_EE+(1i*X1_EE)+1i*Xm_EE*((R2p_EE./st_max_EE)+1i*X2p_EE)./((R2p_EE./st_max_EE)+1i*(Xm_EE+X2p_EE));
I_tmax_EE = V1_EE./Z1_tmax_EE;

Z1_tmax_SE = R1_SE+(1i*X1_SE)+1i*Xm_SE*((R2p_SE./st_max_SE)+1i*X2p_SE)./((R2p_SE./st_max_SE)+1i*(Xm_SE+X2p_SE));
I_tmax_SE = V1_SE./Z1_tmax_SE;


fprintf('\n3.c) Stator current at maximum torque\n');
fprintf('Stator current for SE Motor at maximum torque: %.4f < %.2f° A\n', abs(I_tmax_EE), rad2deg(angle(I_tmax_EE)));
fprintf('Stator current for SE Motor at maximum torque: %.4f < %.2f° A\n', abs(I_tmax_SE), rad2deg(angle(I_tmax_SE)));
disp('At maximum torque, the stator current is higher as the motor draws more current to maintain the required torque. The value depends on the supply voltage and the stator impedance.');

% d) Stator current under no-load conditions (s=0)
Z1_nl_EE = R1_EE+1i*(X1_EE+Xm_EE);
I_nl_EE = V1_EE./Z1_nl_EE;

Z1_nl_SE = R1_SE+1i*(X1_SE+Xm_SE);
I_nl_SE = V1_SE./Z1_nl_SE;

fprintf('\n3.d) Stator current under no-load conditions\n');
fprintf('Stator current for SE Motor under no-load:: %.4f < %.2f° A\n', abs(I_nl_EE), rad2deg(angle(I_nl_EE)));
fprintf('Stator current for SE Motor under no-load:: %.4f < %.2f° A\n', abs(I_nl_SE), rad2deg(angle(I_nl_SE)));
disp('Under no-load conditions, the stator current is relatively low as only the magnetizing current is required to maintain the magnetic field.');

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

