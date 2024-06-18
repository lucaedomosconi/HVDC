
tau = 2.0;
Vmax = 1.5e4;
A = 1e-6;
d1 = 0.5e-4;
d2 = 0.5e-4;

sigma_pap = 1e-14;
sigma_oil = 8.4e-14;
sigma_air = 1e-14;
chi1 = 0.2174; chi2 = 0.5221; chi3 = 3.523; chi4 = 0.4364; chi5 = 1.6996; chi6 = 0.6634;
eps0 = 8.854e-12; eps_pap = 3.14; eps_oil = 2.25;
tau1 = 10; tau2 = 166.7; tau3 = 3333.3; tau4 = 20; tau5 = 600; tau6 = 3333.3;

Qinf = [0,0,0]; Qinfo = Qinf;
Q1 = [0,0,0]; Q2 = [0,0,0]; Q3 = [0,0,0];
Q1o = Q1; Q2o = Q2; Q3o = Q3;
phi = [0];
Cinf = [eps_pap * eps0 * A / d1, eps_oil * eps0 * A / d2];
Ginf = [sigma_pap * A / d1, sigma_oil * A /d2];
C1 = [eps0 * chi1 * A / d1, eps0 * chi4 * A / d2]; C2 = [eps0 * chi2 * A / d1, eps0 * chi5 * A / d2]; C3 = [eps0 * chi3 * A / d1, eps0 * chi6 * A / d2];
%C1 = C1*1e14; C2 = C2*1e14; C3 = C3*1e14;
G1 = [C1(1) / tau1, C1(2) / tau4]; G2 = [C2(1) / tau2, C2(2) /tau5]; G3 = [C3(1) / tau3, C3(2) / tau6];
V = @(t) Vmax * (1 - exp(-t/tau)) * (t <= 30000);

T = 30000;
dt = 0.02;
%time = [0;DATA.time]';
time = 0:dt:T;
I = zeros(1,length(time)-1);
qinfvec = zeros(1,length(time));
q1vec = zeros(1,length(time));
q2vec = zeros(1,length(time));
q3vec = zeros(1,length(time));
%% Homogeneous case
iter = 0;
M = zeros(4,4);
F = zeros(4,1);
M(1,1) = 1;

for t = time(2:end)
    iter = iter+1;
    dt = time(iter+1)-time(iter);

    M(2,2) = 1 + dt / tau1;
    M(3,3) = 1 + dt / tau2;
    M(4,4) = 1 + dt / tau3;

    F(1) = Cinf(1) * V(t);
    F(2) = Q1o(1) + C1(1) * dt / tau1 * V(t);
    F(3) = Q2o(1) + C2(1) * dt / tau2 * V(t);
    F(4) = Q3o(1) + C3(1) * dt / tau3 * V(t);

    sol = M \ F;
    Qinf(1) = sol(1);
    Q1(1) = sol(2);
    Q2(1) = sol(3);
    Q3(1) = sol(4);

    I(iter) = V(t) * Ginf(1) + (Qinf(1) + Q1(1) + Q2(1) + Q3(1) - (Qinfo(1) + Q1o(1) + Q2o(1) + Q3o(1))) / dt;
    Q1o = Q1; Q2o = Q2; Q3o = Q3; Qinfo = Qinf;
end
time = time(2:end);
%% Two layers
iter = 0;
M = zeros(9,9);
F = zeros(9,1);

for t = time(2:end)
    iter = iter+1;
    dt = time(iter+1)-time(iter);
    M(1,1) = 1; M(1,9) = Cinf(1);
    M(2,2) = 1 + dt / tau1; M(2,9) = dt / tau1 * C1(1);
    M(3,3) = 1 + dt / tau2; M(3,9) = dt / tau2 * C2(1);
    M(4,4) = 1 + dt / tau3; M(4,9) = dt / tau3 * C3(1);
    M(5,5) = 1; M(5,9) = - Cinf(2);
    M(6,6) = 1 + dt / tau4; M(6,9) = - dt / tau4 * C1(2);
    M(7,7) = 1 + dt / tau5; M(7,9) = - dt / tau5 * C2(2);
    M(8,8) = 1 + dt / tau6; M(8,9) = - dt / tau6 * C3(2);
    M(9,:) = [-1, -1, -1, -1, 1, 1, 1, 1, dt*(Ginf(2)+Ginf(1))];

    F(1) = Cinf(1) * V(t);
    F(2) = Q1o(1) + C1(1) * dt / tau1 * V(t);
    F(3) = Q2o(1) + C2(1) * dt / tau2 * V(t);
    F(4) = Q3o(1) + C3(1) * dt / tau3 * V(t);
    F(5) = 0;
    F(6) = Q1o(2);
    F(7) = Q2o(2);
    F(8) = Q3o(2);
    F(9) = Qinfo(2) + Q1o(2) + Q2o(2) + Q3o(2) - Qinfo(1) - Q1o(1) - Q2o(1) - Q3o(1) + dt * Ginf(1) * V(t);

    sol = M\F;
    Qinf(1) = sol(1);
    Q1(1) = sol(2);
    Q2(1) = sol(3);
    Q3(1) = sol(4);
    Qinf(2) = sol(5);
    Q1(2) = sol(6);
    Q2(2) = sol(7);
    Q3(2) = sol(8);


    I(iter) = (V(t) - sol(9)) * Ginf(1) + (Qinf(1) + Q1(1) + Q2(1) + Q3(1) - (Qinfo(1) + Q1o(1) + Q2o(1) + Q3o(1))) / dt;
%    I(iter) = sol(9) *Ginf(2);
    Q1o = Q1; Q2o = Q2; Q3o = Q3; Qinfo = Qinf;
end
time = time(2:end);
