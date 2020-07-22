CELLS = table;
[file,path] = uigetfile('*.mat','Choose the combined output files','multiselect','on');
cd(path)
temp = strsplit(path,'/');
Animal = temp{length(temp)-1};
p_struct
information = table; sparsity = table; selectivity = table;
for h = 1:length(file)
    CELLS.path(h,1) = strcat(string(path),file{h});
for g = 1:100
for z = 1:height(CELLS)
    load(CELLS.path(z,1));
    [info,spars,sel] = shuffle_func(Axona_Output,p); 
    infor{1} = info;
    sparsi{1} = spars;
    selec{1} = sel;
end

in = table(infor{1});
sp = table(sparsi{1});
se = table(selec{1});

information = vertcat(information,in);
sparsity = vertcat(sparsity,sp);
selectivity = vertcat(selectivity,se);

process = strcat('Cell_number','_',num2str(z),'_','Shuffle_number','_',num2str(g));
disp(process)
end
end

Shuffle_Results.information = information;
Shuffle_Results.sparsity = sparsity;
Shuffle_Results.selectivity = selectivity;

save('Shuffle_Results.mat','Shuffle_Results');