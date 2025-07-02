figure
hold on
plot(Results.RevRate(isModel), Results.Model_G_a_S_nocut(isModel), 'bo')
plot(Results.RevRate(isModel), Results.Model_G_a_Svd(isModel), 'ro')
plot(Results.RevRate(isModel), Results.Model_G_a_SpSE(isModel), 'go')
legend('no cut)','Vandamme', 'Spherical exp')
hold off
xlabel('RevRate (Myr^{-1})')
ylabel('Model_G_a')