function Default=GetPlotSettingsFun(file)
data=readtable(file);

for m=1:numel(data.Object)
Default.(char(data.Object(m))).Active=data.Active(m);
Default.(char(data.Object(m))).Legend=data.Legend(m);
Default.(char(data.Object(m))).LegendStr=data.LegendStr(m);
Default.(char(data.Object(m))).R=data.B(m);
Default.(char(data.Object(m))).G=data.G(m);
Default.(char(data.Object(m))).B=data.B(m);
Default.(char(data.Object(m))).Color=[data.R(m) data.G(m) data.B(m)];
Default.(char(data.Object(m))).LineWidth=data.LineWidth(m);
Default.(char(data.Object(m))).LineStyle=data.LineStyle(m);
Default.(char(data.Object(m))).FaceAlpha=data.FaceAlpha(m);
end