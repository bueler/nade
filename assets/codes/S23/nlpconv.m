% NLPCONV  Generate and plot convergence for NLP.

mlist = [4, 9, 19, 39, 79, 159];
h = 1 ./ (mlist+1);
fprintf('*** gamma = 0 ***\n')
for j = 1:length(mlist)
    [x,y,U,err0(j)] = nlp(mlist(j),0.0);
end
fprintf('\n*** gamma = 10 ***\n')
for j = 1:length(mlist)
    [x,y,U,err10(j)] = nlp(mlist(j),10.0);
end
loglog(h,err0,'ko'),  hold on
p0 = polyfit(log(h),log(err0),1);
loglog(h,exp(p0(2) + p0(1)*log(h)),'k:')
loglog(h,err10,'k+')
p10 = polyfit(log(h),log(err10),1);
loglog(h,exp(p10(2) + p10(1)*log(h)),'k:')
legend('\gamma=0',sprintf('O(h^{%.2f})',p0(1)),...
       '\gamma=10',sprintf('O(h^{%.2f})',p10(1)),
       'location','northwest')
xlabel h,  ylabel('numerical error')
hold off,  axis tight
