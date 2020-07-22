[file,path] = uigetfile('*.egf','Choose the combined output files','multiselect','on');
cd(path)

if ~iscell(file)
    fi{1} = (file);
    file = fi;
end

for z = 1:length(file)
    [status,eeg_bits{z,1},fs,bytesPerSample] = readEGF(file{z});
end

eeg = vertcat(eeg_bits{:});

fe = strsplit(path,'/');
animal = fe{end-1};
filename = strcat(animal,'_',file{1},'.mat');
cd('/Users/rob/OneDrive - National University of Ireland, Galway/DOWN_SYNDROME_PROJECT');
save(string(filename),'eeg');
cd(path);