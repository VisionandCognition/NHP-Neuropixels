
trials = Log.trial;
for ind = 1:length(trials) %from mat log file
    a = struct2cell(trials(ind))';
    if length(a) < 5
        nC(ind,:) = horzcat(a, 0);
    else
        nC(ind,:) = a;
    end
end

Trials0906R7 = nC;
%Trials1102_r3 = vertcat(nC{:});

jtime = 0.15;
wagjperc = 0.4;

%%noise level 1102
%alpha = 0 
%t = 43:50
%s = 33,34
%alpha = 0.09 
%t = 17,18,23,24,25,26,31,32
%s = 18,20,21,23,26,28,29,31
%alpha = 0.125
%t = 19:22,27:30
%s = 17,19,22,24,25,27,30,32
%alpha = 0.35
%t = 35:42
%s = 1,6,3,8,11,16,9,14

for i=1:length(Trials0906R7)
a=Trials0906R7{i,2};
ind=find(a==33|a==34);
a(ind)=0;
Trials0906R7{i,2}=a;
end

for i=1:length(Trials0906R7)
a=Trials0906R7{i,2};
ind=find(a==18|a==20|a==21|a==23|a==26|a==28|a==29|a==31);
a(ind)=0.15;
Trials0906R7{i,2}=a;
end

for i=1:length(Trials0906R7)
a=Trials0906R7{i,2};
ind=find(a==17|a==19|a==22|a==24|a==25|a==27|a==30|a==32);
a(ind)=0.23;
Trials0906R7{i,2}=a;
end

for i=1:length(Trials0906R7)
a=Trials0906R7{i,2};
ind=find(a==1|a==3|a==6|a==8|a==9|a==11|a==14|a==16);
a(ind)=0.35;
Trials0906R7{i,2}=a;
end
%wrong trials at alpha 0
w0=0; %
for i=1:length(Trials0906R7)
    if isequal(Trials0906R7{i,1},'wrongmatch') & Trials0906R7{i,2}==0
        w0=w0+1;
    end
end

wa0=0;
for i=1:length(Trials0906R7)
    if isequal(Trials0906R7{i,1},'wager') & Trials0906R7{i,2}==0
        wa0=wa0+1;
    end
end

w9=0;
for i=1:length(Trials0906R7)
    if isequal(Trials0906R7{i,1},'wrongmatch') & Trials0906R7{i,2}==0.15
        w9=w9+1;
    end
end

c9=0;
for i=1:length(Trials0906R7)
    if isequal(Trials0906R7{i,1},'correctmatch') & Trials0906R7{i,2}==0.15
        c9=c9+1;
    end
end

wa9=0;
for i=1:length(Trials0906R7)
    if isequal(Trials0906R7{i,1},'wager') & Trials0906R7{i,2}==0.15
        wa9=wa9+1;
    end
end

w125=0;
for i=1:length(Trials0906R7)
    if isequal(Trials0906R7{i,1},'wrongmatch') & Trials0906R7{i,2}==0.23
        w125=w125+1;
    end
end

c125=0;
for i=1:length(Trials0906R7)
    if isequal(Trials0906R7{i,1},'correctmatch') & Trials0906R7{i,2}==0.23
        c125=c125+1;
    end
end

wa125=0;
for i=1:length(Trials0906R7)
    if isequal(Trials0906R7{i,1},'wager') & Trials0906R7{i,2}==0.23
        wa125=wa125+1;
    end
end

w35=0;
for i=1:length(Trials0906R7)
    if isequal(Trials0906R7{i,1},'wrongmatch') & Trials0906R7{i,2}==0.35
        w35=w35+1;
    end
end

c35=0;
for i=1:length(Trials0906R7)
    if isequal(Trials0906R7{i,1},'correctmatch') & Trials0906R7{i,2}==0.35
        c35=c35+1;
    end
end

wa35=0;
for i=1:length(Trials0906R7)
    if isequal(Trials0906R7{i,1},'wager') & Trials0906R7{i,2}==0.35
        wa35=wa35+1;
    end
end

n0 = w0 + wa0;
n9 = c9 + w9 + wa9;
n125 = c125 + w125 + wa125;
n35 = c35 + w35 + wa35;

percc0 = 0;
perf9 = c9/n9;
perf125 = c125/n125;
perf35 = c35/n35;

percwag0 = wa0/n0;
percwag9 = wa9/n9;
percwag125 = wa125/n125;
percwag35 = wa35/n35;

percw0 = w0/n0;
percw9 = w9/n9;
percw125 = w125/n125;
percw35 = w35/n35;

juice0 = jtime*wagjperc*wa0;
juice9 = jtime*c9+jtime*wagjperc*wa9;
juice125 = jtime*c125+jtime*wagjperc*wa125;
juice35 = jtime*c35+jtime*wagjperc*wa35;

figure;
x = categorical({'0','0.15','0.23','0.35'});
y = [percc0,percwag0,percw0;perf9,percwag9,percw9;perf125,percwag125,...
    percw125;perf35,percwag35,percw35];
bar(x,y,'stacked')
ylim([0,1])

% figure;
% x = categorical(["Noise only" "High noise" "Medium noise" "Low noise"]);
% x = reordercats(x,{'Noise only' 'High noise' 'Medium noise' 'Low noise'});
% y = [percwag0,percc0,0,percw0;0,perf9,percwag9,percw9;0,perf125,percwag125,...
%     percw125;0,perf35,percwag35,percw35];
% b = bar(x,y,'stacked');
% %b.FaceColor = 'flat';
% %b(3).CData(1,:) = [0,0.8,0.8];
% ylim([0,1]);
% legend('Correct rejection','Correct match','Opted out','Error');
