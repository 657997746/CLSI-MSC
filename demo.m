clc,clear;
load MSRC_v1.mat  

NMI = [];ACC = [];AR = [];Fscore = [];Precision = [];Recall = [];
for jj = 1:30
    [grps,max_err] = CLSI_MSC(X,gt);
    [A nmi avgent] = compute_nmi(gt,grps);
    acc = Accuracy(grps,double(gt));
    [f,pr,r] = compute_f(gt,grps);
    [ar,RI,MI,HI]=RandIndex(gt,grps);

    NMI = [NMI nmi]; ACC = [ACC acc];
    AR = [AR ar]; 
    Fscore = [Fscore f];
    Precision = [Precision pr];
    Recall = [Recall r];
end
fprintf('NMI of CLSI-MSC is : %.03f(%.03f) ', mean(NMI,2),std(NMI));
fprintf('\nACC of CLSI-MSC is : %.03f(%.03f) ', mean(ACC,2),std(ACC));
fprintf('\nAR of CLSI-MSC is : %.03f(%.03f)  ', mean(AR,2),std(AR));
fprintf('\nFscore of CLSI-MSC is : %.03f(%.03f) ', mean(Fscore,2),std(Fscore));
fprintf('\nPrecision of CLSI-MSC is : %.03f(%.03f) ', mean(Precision,2),std(Precision));
fprintf('\nRecall of CLSI-MSC is : %.03f(%.03f) \n', mean(Recall,2),std(Recall));