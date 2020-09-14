proc contents data=sashelp.bweight;
run;

data bweight;
set sashelp.bweight;
run;

proc print data=bweight (FIRSTOBS = 1 OBS = 20);
run;

proc sgplot data=bweight;
title"Scatterplot (Weight v.s MomWtGain)";
scatter x=MomWtGain y=Weight / markerattrs=(symbol=CircleFilled);
run;

proc corr data=bweight;
var MomWtGain Weight;
run;

proc freq data=bweight;
tables Married*MomSmoke / chisq expected;
run;

title "Vertical Bar Chart";
title2 "Clustered by Married";
proc sgplot data=bweight;
vbar MomSmoke / group=Married groupdisplay=cluster;
run; 



