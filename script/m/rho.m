function rho(level,L)
[t,r]=textread('Arho.txt',"%f %f");
h = L / 2^level;
rho = -r*0.5*h;
figure

plot(t,rho)
xlabel('t [s]')
ylabel('Q [C]')

end
