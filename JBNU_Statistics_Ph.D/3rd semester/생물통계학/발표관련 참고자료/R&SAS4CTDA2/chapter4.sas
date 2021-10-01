/*********************************
	Chapter 4: SAS Programs
********************************/
/*read data from the path of dataset*/
Data dat;
infile "C:/Users/HSY/Desktop/DBP.csv" delimiter="," firstobs=2;
input Subject TRT$ DBP1 DBP2 DBP3 DBP4 DBP5 Age Sex$;
/*create a column diff*/
diff= DBP5-DBP1;
RUN;
proc print data=dat;
run;

/***********************************************
		Section 4.3.1.1
************************************************/
/*test whether the DBP means are different*/
PROC TTEST data= dat;
class TRT;
var DBP1;
RUN;

/*make 2 by 2 table using variable "Sex"*/
PROC FREQ data= dat;
tables Trt*Sex / out= freqs ;
RUN;
proc print data=freqs;
run;

PROC TRANSPOSE data= freqs out= SexbyTrt(drop=_:);
   id  Sex;
   var count;
   by  Trt;
RUN;
/*print the table*/
PROC PRINT data= SexbyTrt;
RUN;

/*test equality of proportions of 2 treatment groups
using Pearson's Chi squares*/
PROC FREQ data= freqs;
weight count;
tables Trt*Sex/ chisq;
RUN;

/*Fit the main effect model on "Sex" and "Age"*/
PROC GLM data= dat;
class Sex(ref="F");
model DBP1= Sex Age / solution;
RUN;

/********************************************
	Section 4.3.1.2
**********************************************/
/*stepwise model selection */
PROC GLMSELECT data= dat;
class Trt Sex;
model diff= Trt|Sex|Age
/ selection= stepwise(select= AIC) stats= all;
RUN;

/*fit the reduced model*/
PROC GLM data= dat;
class TRT(ref="A");
model diff= Age Trt / solution;
RUN;

/******************************
Section 4.3.1.3: MACNOVA for Treatment Difference
*************************************/
/* create the data*/
data mdat;
set dat;
diff2to1 = DBP2-DBP1;
diff3to1 = DBP3-DBP1;
diff4to1 = DBP4-DBP1;
diff5to1 = DBP5-DBP1;
run;
proc print data=mdat;
run;
/* manova using glm*/
PROC glm data= mdat;
class TRT;
model diff2to1 diff3to1 diff4to1 diff5to1= TRT Age/ss3;
contrast '1 vs 2' TRT 1 -1;
manova h=_all_;
RUN;

/*********************************************
		Section 4.3.2
**********************************************/
/*read data from the path of dataset*/
Data dat1;
infile "C:/Users/HSY/Desktop/betablocker.csv" delimiter="," firstobs=2;
input Deaths Total Center Treatment$;
RUN;
proc print data=dat1;
run;

/*fit the logistic regression model*/
PROC GENMOD data= dat1;
class Center(ref="1") Treatment(ref="Control");
model Deaths/Total= Center Treatment/dist= binomial link= logit;
RUN;

/*disperson of parameter*/
/*check table "Criteria For Assessing Goodness Of Fit"
line "pearson chi sqaure*/

/*adjust model fitting with estimate of dispersion*/
PROC GENMOD data= dat1;
class Center(ref="1") Treatment(ref="Control");
model Deaths/Total= Center Treatment
/dist= binomial link= logit scale= pearson;
RUN;

/*fit quasi-likelihood for binomial data*/
PROC GLIMMIX data= dat1;
class Center(ref="1") Treatment(ref="Control");
model Deaths/Total= Center Treatment
		/ link=logit dist=binomial solution;
/*specify the overdispersion parameter*/
random _residual_;
RUN;

/*Section 4.3.3*/
/*read data from the path of dataset*/
Data dat2;
infile "C:/Users/HSY/Desktop/polyps.csv" delimiter="," firstobs=2;
input number treat$ age;
RUN;
proc print data=dat2;
run;

/*fit poisson regression model*/
PROC GENMOD data= dat2;
class treat(ref="placebo");
model number= age|treat / dist= poisson link= log;
RUN;

/*disperson of parameter*/
/*check table "Criteria For Assessing Goodness Of Fit"
line "pearson chi sqaure*/

/*adjust model fitting with estimate of dispersion*/
PROC GENMOD data= dat2;
class treat(ref="placebo");
model number= age|treat / dist= poisson link= log
				scale= PEARSON;
RUN;

/*refit the model without interaction*/
PROC GENMOD data= dat2;
class treat(ref="placebo");
model number= age treat / dist= poisson link= log
				scale= PEARSON;
RUN;

/*disperson of parameter*/
/*check table "Criteria For Assessing Goodness Of Fit"
line "pearson chi sqaure*/

/*fit the quasi poisson*/
PROC GLIMMIX data= dat2;
class treat(ref="placebo");
model number= age treat / dist= poisson link= log solution;
random _residual_;
RUN;

/*fit the negative-binomial model*/
PROC GENMOD data= dat2;
class treat(ref="placebo");
model number= age treat / dist= negbin;
RUN;

/*parameter estimate: k=1/dispersion: 1/0.5816= 1.7194*/
