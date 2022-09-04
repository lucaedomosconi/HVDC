function rho(level,L=0.001)

f = fopen('Arho.txt');
N = textscan(f,'%d',2);
N_v = N{1}(1);
N_t = N{1}(2);

C = textscan(f,'%f');
C = reshape(C{1},[N_v,N_t])';

h = L ./ 2.^level;

figure
xlabel('t [s]')
ylabel('\sigma [C/m^2]')
hold on
grid on
xlim([0,max(C(:,1))])
for i=2:N_v
	plot(C(:,1),C(:,i)*0.5*h(i-1),'DisplayName',strcat('\sigma_',num2str(i-1)))
end
legend('Location','northeast')

end

