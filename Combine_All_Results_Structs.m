[file,path] = uigetfile('*.mat','Choose the combined output files','multiselect','on');
cd(path)
Grand_Results = struct;
for t = 1:size(file,2)
    load(file{t})
    if t == 1
        Grand_Results = Results;
    else
        Grand_Results = vertcat(Grand_Results,Results);
    end
end