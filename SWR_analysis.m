load('/Users/rob/Documents/OneDrive - National University of Ireland, Galway/DOWN_SYNDROME_PROJECT/RawEGFs/ALLEGFS/Ripple_Results.mat');

DS = Results(1:31);
CTL = Results(32:end);

for z = 1:size(DS,2)
    DS_ripples(z,1) = DS(1,z).ripples_per_sec;
    DS_size_ripples(z,1) = size(DS(1,z).ripples,1);
    DS_ripple_freq(z,1) = DS(1,z).ripple_frequency;
    DS_ripple_amp(z,1) = DS(1,z).ripple_amplitude;
end

for z = 1:size(CTL,2)
    CTL_ripples(z,1) = CTL(1,z).ripples_per_sec;
    CTL_size_ripples(z,1) = size(CTL(1,z).ripples,1);
    CTL_ripple_freq(z,1) = CTL(1,z).ripple_frequency;
    CTL_ripple_amp(z,1) = CTL(1,z).ripple_amplitude;
end

DS_mean_ripples_sec = mean(DS_ripples);
CTL_mean_ripples_sec = mean(CTL_ripples);
DS_sem_ripples_sec = sem(DS_ripples);
CTL_sem_ripples_Sec = sem(CTL_ripples);

DS_mean_size_ripples = mean(DS_size_ripples);
CTL_mean_size_ripples = mean(CTL_size_ripples);

for i = 1:size(DS,2)
    temp = struct2table(DS(i).ripples);
    ds_area{i} = (temp.Area);
    clear temp
end

for i = 1:size(CTL,2)
    temp = struct2table(CTL(i).ripples);
    ctl_area{i} = temp.Area;
    clear temp
end


for g = 1:size(ds_area,2)
    if g == 1
        ds_rip_area = table;
        ds_rip_area = ds_area{1};
    else
        ds_rip_area = vertcat(ds_rip_area,ds_area{g});
    end
end

for g = 1:size(ctl_area,2)
    if g == 1
        ctl_rip_area = table;
        ctl_rip_area = ds_area{1};
    else
        ctl_rip_area = vertcat(ctl_rip_area,ds_area{g});
    end
end

for h = 1:size(ds_area,2)
    mean_area_ds(h,1) = mean(ds_area{h});
end

for h = 1:size(ctl_area,2)
    mean_area_ctl(h,1) = mean(ctl_area{h});
end
ds_samples = ds_rip_area/4800;
ctl_samples = ctl_rip_area/4800;
area = vertcat(ds_samples,ctl_samples);
area_grp(1:length(ds_rip_area),1) = {'TS65DN'}; area_grp(length(ds_rip_area)+1:length(ds_rip_area)+length(ctl_rip_area),1) = {'2N Control'};

figure(1)
g = gramm('x',area_grp,'y',area,'group',area_grp,'color',area_grp);
g.geom_jitter();
g.stat_summary('type','sem','geom',{'bar','black_errorbar'});
g.draw();
