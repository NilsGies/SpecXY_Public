function [PeakDB,T]=ReadPeakdata_fun(PeakDBfile)
if istable(PeakDBfile)
    T=PeakDBfile;
else
    T=readtable(PeakDBfile);
end
ColNames=T.Properties.VariableNames;
PeakDB.MineralNames=unique(T.MineralNames);
for n = 1:numel(PeakDB.MineralNames)
    PeakDB.Mineral(n).GroupNames=unique(T.PeakgroupNames(ismember(T.MineralNames,PeakDB.MineralNames(n))));
    for m=1:numel(PeakDB.Mineral(n).GroupNames)
        PeakDB.Mineral(n).Groups(m).PeakNames=unique(T.PeakNames(ismember(T.MineralNames,PeakDB.MineralNames(n))& ismember(T.PeakgroupNames,PeakDB.Mineral(n).GroupNames(m))));
        for l=1:numel(PeakDB.Mineral(n).Groups(m).PeakNames)
            for k=1:numel(ColNames)
                PeakDB.Mineral(n).Groups(m).Peaks(l).(char(ColNames{k}))=T.(char(ColNames{k}))(ismember(T.MineralNames,PeakDB.MineralNames(n))& ismember(T.PeakgroupNames,PeakDB.Mineral(n).GroupNames(m))& ismember(T.PeakNames,PeakDB.Mineral(n).Groups(m).PeakNames(l)));
            end

        end
    end
end
end