%% closed Jung circuit, the effect of arabionose on the green node
clc
clear
close all

%% load data
noara = load('output/leak_0.01_closed_A_100_Ara_0.mat');
ara = load('output/leak_0.01_closed_A_100_Ara_0.001.mat');

%% plot

% initial time for plotting, we will rescale it between 0 and 1
t_init = 100; % min
noara.pos_init = find(noara.t>t_init,1,'first');
ara.pos_init = find(ara.t>t_init,1,'first');

% transformed green concentration
noara.green = (noara.c(noara.pos_init:end,3)-min(noara.c(noara.pos_init:end,3)));
noara.green = noara.green/max(noara.green);
ara.green = (ara.c(ara.pos_init:end,3)-min(ara.c(ara.pos_init:end,3)));
ara.green = ara.green/max(ara.green);

hf = figure;
hold on
plot(noara.t(noara.pos_init:end),noara.green,'LineWidth',2)
plot(ara.t(ara.pos_init:end),ara.green,'LineWidth',2)
xlabel('\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$\hat{c}$ / -','interpreter','latex','Fontsize',18)
% legend('noara','ara','Interpreter','none')
set(gca,'LineWidth',2,'Fontsize',16)
!mkdir -p output
exportgraphics(hf,'output/leak_0.01_closed_A_100_Ara_0.001_compare.pdf','Resolution',300)