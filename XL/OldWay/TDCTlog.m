
Ctrlidx = [Log.trial.Trlnum].'; %might count the last incomplete trial in


% CtrlidxC = find(strcmp({Log.events.type},'correct'));
% Ctrln = [Log.events(CtrlidxC).trialn]'; 

Etrlidx = find(strcmp({Log.events.type},'error'));
Etrln = [Log.events(Etrlidx).trialn]';

Ctrlroll = [1:length(Ctrlidx)-1]'; %correct trials list as extracted by bits

%target distractor location
%1=lower left
%2=upper left
%3=upper right
%4=left right
Targ_idx = [Log.trial.targ_idx].';

LLtarg = find(Targ_idx==1);
ULtarg = find(Targ_idx==2);
URtarg = find(Targ_idx==3);
LRtarg = find(Targ_idx==4);

Dist_idx = [Log.trial.dist_idx].';

LLdist = find(Dist_idx==1);
ULdist = find(Dist_idx==2);
URdist = find(Dist_idx==3);
LRdist = find(Dist_idx==4);

URdist = URdist(1:length(URdist)-1);

LLbg = setdiff(Ctrlroll, union(LLtarg,LLdist));
ULbg = setdiff(Ctrlroll, union(ULtarg,ULdist));
URbg = setdiff(Ctrlroll, union(URtarg,URdist));
LRbg = setdiff(Ctrlroll, union(LRtarg,LRdist));

Motiondir = [Log.trial.motiondir].';

figure;plot(squeeze(nanmean(normMUA(1:50,LLtarg,1:800),[1,2]))); hold on
plot(squeeze(nanmean(normMUA(1:50,LLdist,1:800),[1,2])));
plot(squeeze(nanmean(normMUA(1:50,LLbg,1:800),[1,2])));

figure;
set(gcf,'Position',[0 0 700 500])
h1 = plot(squeeze(nanmean(normMUA(101:150,LLtarg,101:800),[1,2]))); hold on
h2 = plot(squeeze(nanmean(normMUA(101:150,LLdist,101:800),[1,2])));
h3 = plot(squeeze(nanmean(normMUA(101:150,LLbg,101:800),[1,2])));
axis([0,700,-1,1.5])
xticks([0:100:700]);
xticklabels({'-100','0','100','200','300','400','500','600'});
patch([180 390 390 180],[-1 -1 1.5 1.5],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
legend([h1 h2 h3],'Target curve nearby','Distractor curve nearby','Background')
legend('boxoff')


figure;
set(gcf,'Position',[0 0 700 500])
h1 = plot(squeeze(nanmean(normMUA(187:230,LLtarg,101:800),[1,2]))); hold on
h2 = plot(squeeze(nanmean(normMUA(187:230,LLdist,101:800),[1,2])));
h3 = plot(squeeze(nanmean(normMUA(187:230,LLbg,101:800),[1,2])));
legend('Target curve nearby','Distractor curve nearby','Background')
axis([0,700,-1,1.5])
xticks([0:100:700]);
xticklabels({'-100','0','100','200','300','400','500','600'});
patch([180 390 390 180],[-1 -1 1.5 1.5],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')
legend([h1 h2 h3],'Target curve nearby','Distractor curve nearby','Background')
legend('boxoff')



figure;plot(squeeze(nanmean(normMUA(201:230,LLtarg,1:800),[1,2]))); hold on
plot(squeeze(nanmean(normMUA(201:230,LLdist,1:800),[1,2])));
plot(squeeze(nanmean(normMUA(201:230,LLbg,1:800),[1,2])));

figure;plot(squeeze(nanmean(normMUA(300:330,LLtarg,1:800),[1,2]))); hold on
plot(squeeze(nanmean(normMUA(300:330,LLdist,1:800),[1,2])));
plot(squeeze(nanmean(normMUA(300:330,LLbg,1:800),[1,2])));

TARGmean1 = squeeze(nanmean(normMUA(101:150,LLtarg,101:800),[1,2]));
DISTmean1 = squeeze(nanmean(normMUA(101:150,LLdist,101:800),[1,2]));
BGmean1 = squeeze(nanmean(normMUA(101:150,LLbg,101:800),[1,2]));

TARGmDIST1 = TARGmean1 - DISTmean1;
figure;plot(TARGmDIST1)
axis([0,700,-0.15,0.3])
xticks([0:100:700]);
xticklabels({'-100','0','100','200','300','400','500','600'});
patch([180 390 390 180],[-1 -1 1.5 1.5],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')

DISTmBG1 = DISTmean1 - BGmean1;
figure;plot(DISTmBG1)
axis([0,700,-0.15,0.3])
xticks([0:100:700]);
xticklabels({'-100','0','100','200','300','400','500','600'});
patch([180 390 390 180],[-1 -1 1.5 1.5],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')

TARGmean2 = squeeze(nanmean(normMUA(187:230,LLtarg,101:800),[1,2]));
DISTmean2 = squeeze(nanmean(normMUA(187:230,LLdist,101:800),[1,2]));
BGmean2 = squeeze(nanmean(normMUA(187:230,LLbg,101:800),[1,2]));

TARGmDIST2 = TARGmean2 - DISTmean2;
figure;plot(TARGmDIST2)
axis([0,700,-0.3,0.3])
xticks([0:100:700]);
xticklabels({'-100','0','100','200','300','400','500','600'});
patch([180 390 390 180],[-1 -1 1.5 1.5],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')

DISTmBG2 = DISTmean2 - BGmean2;
figure;plot(DISTmBG2)
axis([0,700,-0.3,0.3])
xticks([0:100:700]);
xticklabels({'-100','0','100','200','300','400','500','600'});
patch([180 390 390 180],[-1 -1 1.5 1.5],[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeColor','none')

figure;
for i=101:149
    subplot(7,7,i-100)
    plot(squeeze(nanmean(normMUA(i,LLtarg,1:800)))); hold on
    plot(squeeze(nanmean(normMUA(i,LLdist,1:800))));
    plot(squeeze(nanmean(normMUA(i,LLbg,1:800))));
end

figure;
for i=201:249
    subplot(7,7,i-200)
    plot(squeeze(nanmean(normMUA(i,LLtarg,1:800)))); hold on
    plot(squeeze(nanmean(normMUA(i,LLdist,1:800))));
    plot(squeeze(nanmean(normMUA(i,LLbg,1:800))));
end

figure;
for i=1:100
    subplot(10,10,i)
    plot(squeeze(mean(normMUA(i,LLtarg,101:900))))
end

figure;
for i=1:100
    subplot(10,10,i)
    plot(squeeze(mean(normMUA(i,LRtrial,:))))
end
% 
figure;
for i=1:length(LLtarg(1:100))
    subplot(10,10,i)
    plot(squeeze(normMUA(111,LLtarg(i+80),1:1000)))
end
figure;
for i=1:length(LLtarg(1:100))
    subplot(10,10,i)
    plot(squeeze(normMUA(306,LLtarg(i+80),1:800)))
end