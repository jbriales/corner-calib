function plotLIDARframe()
xlabel('X [m]')
ylabel('Y [m]')
plot(0,0,'>', 'MarkerSize',10, 'LineWidth',3)
plotHomLineWin( [1 0 0], 'k' );
plotHomLineWin( [0 1 0], 'k' );
end
