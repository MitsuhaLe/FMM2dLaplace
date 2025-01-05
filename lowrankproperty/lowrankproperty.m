%need smatrix function

ntarget = 500;
x = rand(2, ntarget);
nsource = 1000;
y = rand(2, nsource);
for gap = 0.1:0.05:2
    yGap = y + [gap+1; 0];
    S = smatrix(x, yGap);
    plot(gap, rank(S),'k.');
    hold on
end
ylabel("rank")
xlabel("gap")
hold off

figure(2)
plot(x(1, :), x(2, :),'b.');
hold on
plot(y(1, :) + 1.5, y(2, :), 'k.');
axis equal
text(0.5, -0.2, "target:X","FontSize",14);
text(2, -0.2, "source:Y","FontSize",14);
hold off