% EIGPATH  Solve P23 a) on Assignment #5.

M = @(x) [2,x,x; -1,0,1; -1,1,0];
figure(1),  hold on
for x=linspace(-1,5,301)
    lam = eig(M(x));
    [lamx,yi] = sort(real(lam));
    lamy = imag(lam(yi));
    if floor(x) == x
        MS = 18.0;
    else
        MS = 10.0;
    end
    plot(lamx(1),lamy(1),'b.','markersize',MS)
    plot(lamx(2),lamy(2),'g.','markersize',MS)
    plot(lamx(3),lamy(3),'r.','markersize',MS)
    if floor(x) == x
        for j=1:3
            text(lamx(j)+0.1,lamy(j)+0.2,num2str(x),...
                 'color','k','fontsize',14.0);
        end
    end
end
hold off
axis([-2 4 -5 5]),  grid on
title('eigenvalues of M(x) for different x values')
