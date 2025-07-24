
muafile = '\\vs03\VS03-VandC-6\Neuropixels_NHP\Data_collection\m032\20240828\recording1\tempmua';

normname = '\\vs03\VS03-VandC-6\Neuropixels_NHP\Data_collection\m032\20240828\recording1\normmua';

load(muafile);

trial_length = 0.6;
pre_trial = 0.2;
post_trial = trial_length-pre_trial;
tb =-pre_trial*1000:1:(post_trial)*1000;

baseT = tb > -100 & tb <= 0;

normMUA = nan(size(temp_mua));
for chn = 1:384
    temp_chn_all = [];
    base = nanmean(nanmean(temp_mua(chn,:,baseT),3),2);
    temp_trial_mua = squeeze(temp_mua(chn,:,:))-base;
    noise = nanstd(nanmean(temp_trial_mua(:,baseT),2));
    signal = smooth(squeeze(nanmean(temp_trial_mua(:,tb>0.02))),25,'lowess');
    temp_trial_mua = temp_trial_mua./max(signal);
    temp_chn_all = [temp_chn_all; temp_trial_mua];
    maxsignal = max(signal);
    SNR(chn) = maxsignal/noise;
    normMUA(chn,:,:) = temp_chn_all;
    chn;
end

save(normname,'normMUA','SNR','tb','-v7.3')