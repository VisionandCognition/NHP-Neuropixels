
%import stimulus parameters from logfile
logdir = '\\vs03\VS03-VandC-6\Neuropixels_NHP\Data_collection\m032\20240904\recording2\run-002_T-151500\';
logfile = dir([logdir,'*.mat']).name;
savename = '\\vs03\VS03-VandC-6\Neuropixels_NHP\Data_collection\m032\20240904\recording2\barmap';

load([logdir,logfile]);

saveout = 1;

analyzedata = 0;

% To plot dots on top of data:
%Useful for testing potential positions for stimuli
%Set to empty e.g. RFx = []; if you don;t want this
RFx = []; %
RFy = []; % 

fig = 1;

notch = 1;

%split channels into sensible banks, or not
bank = zeros(1,384);
bank(1:48) = 1;
bank(49:96) = 2;
bank(97:144) = 3;
bank(145:192) = 4;
bank(193:240) = 5;
bank(241:288) = 6;
bank(289:336) = 7;
bank(337:384) = 8;
bankcol = [1 0 0;0 1 0;0 0 1;0 0 0;1 1 0;0 1 1;1 0 1;1 0.5 0];
% bank(1:150) = 1;
% bank(151:384) = 2;
% bankcol = [1 0 0;0 1 0];

speed = Stm.BarSpeed;
bardur = Stm.BarDur;
bardist = Stm.BarDist;
xo = Stm.RFx;
yo = Stm.RFy;
pixperdeg = Par.PixPerDeg;

trlcidx = find(strcmp({Log.events.type},'correct'));
trldir = [Log.events(trlcidx).dir];

%0 = left to right, 90 = down-to-up, 180 = right to left, 270 = up to down
trl1 = find(trldir==0);
trl2 = find(trldir==90);
trl3 = find(trldir==180);
trl4 = find(trldir==270);
trlid = trldir;
trlid(trl1) = 1;
trlid(trl2) = 2;
trlid(trl3) = 3;
trlid(trl4) = 4; %assign trial id

%adjust based on actual trial information
Fs = 30000;
trial_length = 2;
pre_trial = 0.2;
post_trial = trial_length-pre_trial;
tb =-pre_trial*1000:1:(post_trial)*1000;


nchans = size(temp_mua,1);
ModS = zeros(nchans,4);
ModM = zeros(nchans,4);
Ons = zeros(nchans,4);
Offs = zeros(nchans,4);

clear SNR
if analyzedata
    if notch
        mon = 60; %Frequency of minitor
        no = 2;
        wn = [mon-2 mon+2]./(Fs./2);
        [b1,a1] = butter(no,wn,'stop');
        wn = [mon*2-2 mon*2+2]./(Fs./2);
        [b2,a2] = butter(no,wn,'stop');
        wn = [mon*3-2 mon*3+2]./(Fs./2);
        [b3,a3] = butter(no,wn,'stop');
    end

    if fig
        figure,hold on
        colind = jet(96);
    end

    for ch = 1:size(temp_mua,1)
        clear mua
        mua = [];
        buf = squeeze(temp_mua(ch,:,:));

        if notch
            %Filter out the monitor frequency
            buf = filtfilt(b1,a1,buf);
            buf = filtfilt(b2,a2,buf);
            buf = filtfilt(b3,a3,buf);
        end
        mua = [mua,buf];

        for n = 1:4
            %Get trials with this motion direction
            f = find(trlid == n);

            %Average them
            MUAm(n,:) = nanmean(mua(f,:));
            MUAs(n,:) = nanstd(mua(f,:))./sqrt(length(f));

            %Get noise levels before smoothing
            BaseT = find(tb >-100 & tb < 0);
            Base = nanmean(MUAm(n,BaseT));
            BaseS = nanstd(MUAm(n,BaseT));

            %Smooth it to get a maximum...
            gt = find(tb>0 & tb<2000);
            sm = smooth(MUAm(n,gt),10);
            [mx,mi] = max(sm);
            Scale = mx-Base;

            %Is the max significantly different to the base?
            SigDif(ch,n) = mx > (Base+(1.*BaseS));


             %Automatic RF fitting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %Fit a guassian to it...
             %Get fitting ranges
             RespT = find(tb > 0);
             %time
             td = tb(RespT);
             index = 1:5:length(td);
             pxm = td(index);
             sdind = [10:10:200];

             z = 0;
             for i = 1:length(pxm)
                 for j = 1:length(sdind)
                     norm = normpdf(td,pxm(i),sdind(j));
                     norm = norm./(max(norm));
                     model = (norm.*Scale)+Base;
                     diff = sum((model-MUAm(n,RespT)).^2);
                     z = z+1;
                     ind(z,:) = [diff,pxm(i),sdind(j)];
                 end
             end
             [x,y] = min(ind(:,1));
             ModM(ch,n) = ind(y,2);
             ModS(ch,n) = ind(y,3);
             Ons(ch,n) = ModM(ch,n)-(1.95.*ModS(ch,n));
             Offs(ch,n) = ModM(ch,n)+(1.95.*ModS(ch,n));
             %Bet model
             norm = normpdf(td, ModM(ch,n), ModS(ch,n));
             norm = norm./(max(norm));
             bestmodel = (norm.*Scale)+Base;


            if fig
                subplot(2,2,n),plot(td,MUAm(n,RespT),'b',td,bestmodel,'r')
            end


            if fig
                h = line([Ons(ch,n),Ons(ch,n)],get(gca,'YLim'));
                set(h,'Color',[1 0 1])
                h = line([Offs(ch,n),Offs(ch,n)],get(gca,'YLim'));
                set(h,'Color',[1 0 1])
            end

            SNR(ch,n)=Scale/BaseS;
        end
        ch

    end

    if saveout
    save(savename,'Ons','Offs','SigDif','speed','SNR','Log','pixperdeg')
    end
else
    load(savename)
end
    
    
%SKIP TO HERE
SNR = mean(SNR,2);
SNRcutoff = 13;
figure;
set(gcf,'Position',[0 0 1024 768])
for ch = 1:size(temp_mua,1)
    %ONly plot channels where all directon were signifcant and teh SNR is
    %high enough
    if sum(SigDif(ch,:)) == 4 && SNR(ch)>SNRcutoff %SNR(ch)>SNRcutoff
        %Now distance = speed*time
        %This gives distanbce travelled by bar in pixels before the onset and
        %offset
        onsdist = speed.*(Ons(ch,:)./1000);
        offsdist = speed.*(Offs(ch,:)./1000);
        
        %Stimuli 1-4 go
        %1 = horizontal left-to-right (180 deg),
        %2 = bottom to top 270
        %3 = right to left 0
        %4 = top to bottom 90
        angles = [180 270 0 90];
        
        %Get starting position of bars
        sx = xo+(bardist./2).*cosd(angles);
        sy = yo+(bardist./2).*sind(angles);
        
        %Angular distance moved
        %(direction is opposite to angle of starting
        %position)
        on_angx = onsdist.*cosd(180-angles);
        on_angy = onsdist.*sind(angles);
        off_angx = offsdist.*cosd(180-angles);
        off_angy = offsdist.*sind(angles);
        
        %So the on and off points are starting position + angular distance...
        onx = sx+on_angx;
        ony = sy-on_angy;
        offx = sx+off_angx;
        offy = sy-off_angy;
        
        %get RF vboundaries
        bottom = (ony(2)+offy(4))./2;
        right = (onx(1)+offx(3))./2;
        top =   (ony(4)+offy(2))./2;
        left =   (onx(3)+offx(1))./2;
        
        RF.centerx(ch) = (right+left)./2;
        RF.centery(ch) = (top+bottom)./2;
        
        RF.sz(ch) = sqrt(abs(top-bottom).*abs(right-left));
        RF.szdeg(ch) = sqrt(abs(top-bottom).*abs(right-left))./pixperdeg;
        
        XVEC1 = [left  right  right  left  left];
        YVEC1 = [bottom bottom  top top  bottom];
        
        RF.XVEC1(ch,:) = XVEC1;
        RF.YVEC1(ch,:) = YVEC1;
        
        h = line(XVEC1,YVEC1);
        set(h,'Color',bankcol(bank(ch),:))
%         axis square
        hold on
        scatter(0,0,'r','f')
        scatter(RF.centerx(ch),RF.centery(ch),36,bankcol(bank(ch),:),'f')
        axis([-512 512 -384 384])
        xticks([-20*pixperdeg -15*pixperdeg -10*pixperdeg -5*pixperdeg 0 ...
            5*pixperdeg 10*pixperdeg 15*pixperdeg 20*pixperdeg])
        xticklabels({'-20','-15','-10','-5','0','5','10','15','20'})
        yticks([-15*pixperdeg -10*pixperdeg -5*pixperdeg 0 ...
            5*pixperdeg 10*pixperdeg 15*pixperdeg])
        yticklabels({'-15','-10','-5','0','5','10','15'})
        hold on,scatter(sx,sy)
%         disp(['channel: ' ,num2str(ch)])
%         disp(['centerx = ',num2str(RF.centerx(ch))])
%         disp(['centery = ',num2str(RF.centery(ch))])
        %position in degrees
        RF.ang(ch)= atand(RF.centery(ch)./RF.centerx(ch));
        
        %pix2deg conversion
        RF.ecc(ch) = sqrt(RF.centerx(ch).^2+RF.centery(ch).^2)./pixperdeg;
        
        % disp(['Angle = ',num2str(RF_ang(ch))])
        disp(['Ecc = ',num2str(RF.ecc(ch))])
        disp(' ')
        
        text(RF.centerx(ch),RF.centery(ch),num2str(ch))
        %Save out centx
    else
        %Bad channels get set to NaN
        RF.centerx(ch)=NaN;
        RF.centery(ch)=NaN;
        RF.sz(ch)=NaN;
        RF.szdeg(ch)=NaN;
        RF.XVEC1(ch,:)=NaN(1,5);
        RF.YVEC1(ch,:)=NaN(1,5);
        RF.ang(ch) = NaN;
        RF.ecc(ch) = NaN;
    end
end

if ~isempty(RFx)
    %SCatter on markers
    hold on,scatter(RFx,RFy,'MarkerFaceColor',[0.8 0.8 0.8])
    for i = 1:length(RFx)
        text(RFx(i),RFy(i),(['x=' num2str(RFx(i)) ', y=' num2str(RFy(i))]))
    end
end

figure;
for a = 1:8
    subplot(2,4,a),hold on
    for ch = find(bank==a)
        h = line(RF.XVEC1(ch,:)./pixperdeg,RF.YVEC1(ch,:)./pixperdeg);
        set(h,'Color',bankcol(a,:))
        scatter(RF.centerx(ch)./pixperdeg,RF.centery(ch)./pixperdeg,36,bankcol(a,:),'f')
        text(RF.centerx(ch)./pixperdeg,RF.centery(ch)./pixperdeg,num2str(ch))
    end
    scatter(0,0,'r','f')
    axis square
    axis([-20 20 -20 20])
end
