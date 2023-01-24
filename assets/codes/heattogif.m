% HEATTOGIF  Generate animated .gif from HEAT output.
% author: Stefano F.
% requires: Matlab and HEAT
% note: For generating .gif in Octave see
%       https://electroagenda.com/en/create-gif-files-in-octave-and-matlab/)

[t, x, T] = heat(20,20);
[m, n] = size(T);
xx = x(1,:);
p = plot(xx,T(1,:));
axis([0 1 0 1])
for i = 1:m
    p.YData = T(i,:);
    exportgraphics(gca,"heat.gif","Append",true)
end
