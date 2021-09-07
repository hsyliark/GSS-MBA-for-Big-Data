/************************************
  Section 1.4.1.2: Data Generation
*************************************/
/* simulate the "Placebo"  */
data dat4placebo(keep =age bp4base bp4end bp4diff trt);
trt = "Placebo";
n=100;mu=100;sd=10;mud =20;agemu=50;agesd=10;
do i =1 to n;
age     = agemu+agesd*rand("Normal");
bp4base = mu+sd*rand("Normal");
bp4end  = mu+sd*rand("Normal");
bp4diff = bp4end-bp4base;
output;
end;
run;
/* simulate the "drug"  */
data dat4drug(keep =age bp4base bp4end bp4diff trt);
trt = "Drug";
n=100;mu=100;sd=10;mud =20;agemu=50;agesd=10;
do i =1 to n;
age     = agemu+agesd*rand("Normal");
bp4base = mu-mud+sd*rand("Normal");
bp4end  = mu+sd*rand("Normal");
bp4diff = bp4end-bp4base;
output;
end;
run;
/* combine these two datasets into "dat"*/
data dat;
set dat4placebo dat4drug;
run;

proc print data=dat;
run;

/*********************************************
  Section 1.4.2 Data analysis
**********************************************/
ods graphics on;
proc glm data=dat;
class trt(ref="Placebo");
model bp4diff = age trt trt*age/solution;
run;
ods graphics off;


/*Section 3.3.1.1*/
/*read data*/
Data dat;
/*read data from the path of dataset*/
infile "Your file path/DBP.csv" delimiter="," firstobs=2;
input Subject TRT$ DBP1 DBP2 DBP3 DBP4 DBP5 Age Sex$;
/*create a column "diff"*/
diff= DBP5-DBP1;
RUN;

/*show the first 6 observations using proc print*/
PROC PRINT data=dat(obs=6);
RUN;

/*Section 3.3.1.2*/
/*call two-side t-test with equal variance*/
PROC TTEST data=dat sides=2 alpha=0.05 h0=0;
/*specifies the title line of result*/
title "Two sample two side t-test";
/*defines the grouping variable*/
class TRT;
/*variable whose means will be compared*/
var diff;
RUN;

/*welch t-test with unequal variance*/
/*check the lines labeled "Satterthwaite" from above results*/

/*test the null hypothesis for equal variance*/
/*check the table named "Equality of Variance" from above results*/

/*nonparametric t-test(Wilcoxon rank-sum test)*/
PROC NPAR1WAY wilcoxon data=dat;
title "Wilcoxon rank-sum test";
class TRT;
var diff;
exact wilcoxon;                 /*request for exact p-value*/
RUN;

/*call one-side t-test */
/*alternative hypothesis: means(A)<means(B)*/
PROC TTEST data=dat sides=L alpha=0.05 h0=0;
title "Two sample one-side t-test";
class TRT;
var diff;
RUN;

/*Section 3.3.1.3*/
/*bootstrap method*/
/*call bootstrap with 1000 replication*/
PROC SURVEYSELECT data=dat out=boot_samples seed=123
/* specify the type of random sampling */
method=urs
/* get a sample of the same size as original data set */
samprate=1
/*give the times a record chosen*/
outhits             /
* specify the number of bootstrap samples */
rep=1000;
/*bootstrapping by TRT without overlap*/
strata TRT / alloc = proportional;
RUN;

/*calculate the means of each replication*/
PROC MEANS data=boot_samples mean;
var diff;
class Replicate TRT;
/*output the results*/
output out=temp1(where=(_type_=3) drop= _freq_) Mean=mean;
RUN;

/*transform the data to contingency table*/
PROC TRANSPOSE data=temp1(drop=_type_) out=temp2(drop=_:);
id Trt;
var Mean;
by Replicate;
RUN;

/*calculate the mean difference*/
Data boot_diff;
set temp2;
mean_diff=A-B;
drop A B;
RUN;

/*calculate and show the bootstraps quantiles*/
PROC UNIVARIATE data= boot_diff;
var mean_diff;
output out= result PCTLPTS= 2.5 50 97.5 PCTLPRE = q;
RUN;

/*Section 3.3.1.4*/
/*extract the means by treatment group*/
PROC MEANS data= dat;
/*defines the grouping variable*/
class TRT;
/*variable whose means will be extracted*/
var DBP1 DBP2 DBP3 DBP4 DBP5;
RUN;

/*rearrange the dataframe of dat into a "long" format*/
DATA Dat_L;
set dat;
/*create an array*/
array ADBP(1:5) DBP1-DBP5;
/*set new variable (DBP) equal to the value of
                 the array for the given time*/
do time= 1 to 5;
DBP= ADBP(time);
/*output the results*/
output;
end;
/*drop unused variables*/
drop DBP1-DBP5 diff;
RUN;

/* sort data in ascending order by time*/
PROC SORT data= Dat_L out= Dat_L;
by time;
RUN;

/*show the first 6 observations*/
PROC PRINT data= Dat_L(obs=6);
RUN;

/*test treatment "A"*/
DATA datA(where=(TRT='A'));
set Dat_L;
RUN;

PROC ANOVA data= datA;
title "one-way ANOVA to test treatment 'A'";
class time;
model DBP= time;
/*determines the times at which DBP
     means are significantly different*/
means time/ tukey cldiff;
RUN;

/*test treatment "B"*/
DATA datB(where=(TRT='B'));
set Dat_L;
RUN;

PROC ANOVA data= datB;
title "one-way ANOVA to test treatment 'B'";
class time;
model DBP= time;
means time / tukey cldiff;
RUN;

/*determines the times at which DBP means are significantly different*/
/*check the result from table named
"Tukey's Studentized Range (HSD) Test for DBP"*/

/*Section 3.3.1.5*/
/*two-way anova test the significance of interaction using glm*/
PROC GLM data= Dat_L;
title "two-way ANOVA test using glm";
class TRT time;
model DBP= TRT time TRT*time;
lsmeans TRT time TRT*time/ pdiff adjust= tukey;
RUN;

/*Section 3.3.2.1*/
DATA Ulcer;                   /*load data*/
input Trt$ x_4 n @@;
Heal="Yes"; Count=x_4; output;
Heal="No"; Count=n-x_4; output;
datalines;
d1 69 168 d2 113 182 d3 120 165 d4 145 188
RUN;

/*test equal probabilities of healing for the four treatment groups
using Pearson's chi square test*/
PROC FREQ data= Ulcer order=data;
weight Count;
tables Trt*Heal/ chisq;
RUN;

/*comparison between two treatment groups*/
/*load data of group1(0 mg)and group3(800 mg)*/
DATA dat_1_1(where=(Trt="d1"|Trt="d3"));
Set Ulcer;
RUN;
PROC FREQ data= dat_1_1 order=data;
weight Count;
tables Trt*Heal/ chisq;
RUN;

/*load data of group2(400 mg)and group3(800 mg)*/
DATA dat_1_2(where=(Trt="d2"|Trt="d3"));
Set Ulcer;
RUN;
PROC FREQ data= dat_1_2 order=data;
weight Count;
tables Trt*Heal/ chisq;
RUN;

/*load data of group3(800 mg)and group4(1600 mg)*/
DATA dat_1_3(where=(Trt="d3"|Trt="d4"));
Set Ulcer;
RUN;
PROC FREQ data= dat_1_3 order=data;
weight Count;
tables Trt*Heal/ chisq;
RUN;


/*transform dataframe to contingency table*/
/*sort data in ascending order by Trt*/
PROC SORT data= Ulcer out=tab_Ulcer;
by Trt;
RUN;

PROC TRANSPOSE data=tab_Ulcer out=tab_Ulcer(drop=_:);
   id  Heal;
   var Count;
   by  Trt;
RUN;


/*********************************
	Chapter 4: SAS Programs
********************************/
/*read data from the path of dataset*/
Data dat;
infile "Your data path/DBP.csv" delimiter="," firstobs=2;
input Subject TRT$ DBP1 DBP2 DBP3 DBP4 DBP5 Age Sex$;
/*create a column diff*/
diff= DBP5-DBP1;
RUN;

/***********************************************
		Section 4.3.1.1
************************************************/
/*test whether the DBP means are different*/
PROC TTEST data= dat;
title "Two sample two sided t-test";
class Trt;
var DBP1;
RUN;

/*make 2 by 2 table using variable "Sex"*/
PROC FREQ data= dat;
tables Trt*Sex / out= freqs ;
RUN;

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
infile "your path/betablocker.csv" delimiter="," firstobs=2;
input Deaths Total Center Treatment$;
RUN;

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
infile "your path/polyps.csv" delimiter="," firstobs=2;
input number treat$ age;
RUN;

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


/***************************************************
	Chapter 5
********************************************************/
/*read data from the path of dataset*/
Data dat3;
infile "your path/CTCarcinoma.csv" delimiter="," firstobs=2;
input TRT$ Time Status Age;
RUN;

/*print the first 3 observations*/
PROC PRINT data= dat3(obs=3);
RUN;

PROC SORT data= dat3 out= dat3;
by descending TRT;
RUN;

/******************************************************
	Section 5.5.1.1
*********************************************************/
/*fit Kaplan-Meier and
estimat suvival and cumulative harzard function and plot*/
PROC LIFETEST data= dat3;
time Time*Status(1);       /*0 is event, 1 is censored*/
strata TRT;
RUN;

/*test the treatment difference*/
/*check table "Stratified Test of Equality over Group"*/

/*******************************************************
	Section 5.5.1.2
***********************************************************/
/*fit exponential model*/
PROC LIFEREG data= dat3 order= data;
class TRT;
model Time*Status(1)= TRT / dist= exponential;
RUN;

/*fit Weibull model*/
PROC LIFEREG data= dat3 order= data;
class TRT;
model Time*Status(1)= TRT / dist=weibull;
RUN;

/*fit exponential model+Age*/
PROC LIFEREG data= dat3 order= data;
class TRT;
model Time*Status(1)= TRT Age / dist= exponential;
RUN;

/*fit Weibull model+Age*/
PROC LIFEREG data= dat3 order= data;
class TRT;
model Time*Status(1)= TRT Age / dist= weibull;
RUN;

/*************************************************
	Section 5.5.1.3
*************************************************/
/*fit Cox regression model*/
PROC PHREG data= dat3;
class TRT(ref= "S+CT");
model Time*Status(1)= TRT/ ties=efron;
RUN;

/*fit Cox regression model+Age*/
PROC PHREG data= dat3;
class TRT(ref= "S+CT");
model Time*Status(1)= TRT Age/ ties=efron;
RUN;

/************************************************
	Section 5.5.2
************************************************/
/*read data from the path of dataset*/
Data dat4;
infile "/home/ln23040/BreastCancer.csv" delimiter="," firstobs=2;
input tL tU$ TRT Status;
if tU= "NA" then tU= .;
ntU= tU+0;      /* transform the type of variable "tU"
                    from char to int*/
if ntU= . then time= tL;
else time= (tL+ntU)/2;
if tL= 0 then ntL= .;
else ntL= tL;
RUN;

PROC SORT data= dat4 out= dat4;
by descending TRT;
RUN;

/*************************************************
	Section 5.5.2.1
**************************************************/
/*fit Turnbull?s estimator*/
PROC ICLIFETEST data= dat4 impute(seed= 123) method= turnbull;
strata TRT;
time (tL,ntU);
RUN;

/*fit Kaplan-Meier estimator with midpoint or left point*/
PROC LIFETEST data= dat4;
time time*Status(0);         /*0 is censored, 1 is event*/
strata TRT;
RUN;

/*****************************************************
	Section 5.5.2.2
********************************************************/
/*fit exponential model*/
PROC LIFEREG data= dat4 order= data;
class TRT;
model (ntL, ntU)= TRT / dist= exponential;
RUN;

/*fit Weibull model*/
PROC LIFEREG data= dat4 order= data;
class TRT;
model (ntL, ntU)= TRT / dist= weibull;
RUN;

/*************************************
	Section 6.3.1.1
*****************************************/
/*read data from the path of dataset*/
Data dat;
infile "your path/DBP.csv" delimiter="," firstobs=2; 
input Subject TRT$ DBP1 DBP2 DBP3 DBP4 DBP5 Age Sex$;
diff= DBP5-DBP1;       /*create a column "diff"*/
RUN;

/*rearrange the dataframe of dat into a "long" format*/
DATA Dat_L;
set dat;
array ADBP(1:5) DBP1-DBP5;  /*creat an array*/
/*set new variable (DBP) equal to the value of 
               the array for the given time*/
do time= 1 to 5;  
DBP= ADBP(time);
output;                 /*output the results*/
end;
drop DBP1-DBP5 diff;    /*drop unuseful variables*/
RUN;

PROC SORT data= Dat_L out=Dat_L1; 
by time;
RUN;

/*show the first 6 observations*/
PROC PRINT data= Dat_L1(obs=6);
RUN;

/*fit linear regression model for 40 patients and 
                        extract intercept and slope*/
PROC REG data= Dat_L outest= coef;
by Subject;
model DBP= Time;
RUN;

/*print out the coefficient:
  subject #1-20 are A and the others are B*/
PROC PRINT data= coef;     
RUN;

/*add column "Trt" into table coefficient*/
PROC SQL;
CREATE TABLE dat_coef AS
SELECT dat.Subject, dat.TRT, coef.Intercept, coef.time
FROM dat FULL JOIN coef
ON dat.subject= coef.subject;
QUIT;

/*model the slope and intercept relationship by linear regression*/
/*fit model 1 with interaction*/
PROC GLM data= dat_coef;
class TRT(ref= "A");
model time= Intercept TRT Intercept*TRt/ solution;
RUN;

/*fit model 2 without interaction*/
PROC GLM data= dat_coef;
class TRT(ref= "A");
model time= Intercept TRT / solution;
RUN;

/*test difference of intercept and slop between two treatments*/ 
PROC TTEST data= dat_coef sides=2 alpha=0.05 h0=0;
class TRT;                           
var time;                           
RUN;

PROC TTEST data= dat_coef sides=2 alpha=0.05 h0=0;
class TRT;                           
var Intercept;                           
RUN;

/*****************************************************
	Section 6.3.1.2
*********************************************************/
/*fit model 1(set both intercept and slope as random effect)*/
PROC MIXED data=Dat_L  method= ML;
class TRT(ref= "A");
model DBP= TRT Time TRT*Time / s;
random int time / type= un subject=Subject;
RUN;

/*fit model 2(set intercept as random effect)*/
PROC MIXED data= Dat_L method= ML;
class TRT(ref= "A");
model DBP= TRT Time TRT*Time / s;
random int / type= un subject= Subject;
RUN;

/*fit model 3*/
PROC MIXED data= Dat_L;
class TRT(ref= "A");
model DBP= TRT Time / s;
random int time/ type= un subject= Subject;
RUN;

/*fit model 4*/
PROC MIXED data= Dat_L;
class TRT(ref= "A");
model DBP= TRT Time / s;
random int / type= un subject= Subject;
RUN;

/*model 5: fit model 3 include "Age" effect*/
PROC MIXED data= Dat_L method= ML;
class TRT(ref= "A");
model DBP= TRT Time Age / s;
random int time/ type= un subject= Subject;
RUN;


/*model 6: fit model 3 include "Age" and "Sex" effect*/
PROC MIXED data= Dat_L method= ML;
class TRT(ref= "A") Sex(ref= "F");
model DBP= TRT Time Age Sex / s;
random int time/ type= un subject= Subject;
RUN;

/****************************************************
	Section 6.3.2.1
********************************************************/
/*read data from the path of dataset*/
Data dat5;
infile "your path/Ulcer.csv" delimiter="," firstobs=2;  
input Subject TRT WeekH Time0 Time1 Time2 Time4;
RUN;

/*print the first 6 observations*/
PROC PRINT data= dat5(obs= 6);
RUN;

/*summary information for the four treatment groups 
       and four time points*/
PROC FREQ data= dat5;
tables TRT TRT*Time1 TRT*Time2 TRT*Time4;
RUN;

DATA Dat_L2;
set dat5;
array AT[3] Time1 Time2 Time4;   /*create an array*/
do i= 1 to 3;                    /*set new variable (DBP) equal to the value of 
                                   the array for the given time*/
if i= 3 then Time= i+1;
else Time= i;
Heal= AT(i);
output;                          /*output the results*/
end;
drop i Time0 Time1 Time2 Time4;    /*drop non-use variables*/
RUN;


PROC SORT data= Dat_L2 out= Dat_L3;
by Subject;
RUN;

/***************************************************************
	Section 6.3.2.2
****************************************************************/
/*fit model 1: with interaction*/
PROC GENMOD data= Dat_L2 descending;
class Trt(ref="1") Time(ref="1");
model Heal= Trt Time Trt*Time/ dist= binomial link= logit;
RUN;

/*fit model 2: without interaction*/
PROC GENMOD data= Dat_L2 descending;
class Trt(ref="1") Time(ref="1");
model Heal= Trt Time/ dist= binomial link= logit;
RUN;

/*multiple comparison using "Tukey" test*/
PROC GLM data= Dat_L2;
class Trt Time;
model Heal= Trt Time;
means TRT Time / Tukey;
run;

/**************************************************************
	Section 6.3.2.3
**************************************************************/
/*fit model 3: Generalized Linear Mixed Model!!!some problem*/
PROC GLIMMIX data= DAT_L2;
class Trt(ref="1") Time(ref="1");
model Heal= Trt/ dist=binomial link=logit solution;
random int / subject= Subject;
RUN;

/*fit model 4*/
PROC GENMOD data= Dat_L2 descending;
class Trt(ref="1") Time(ref="1");
model Heal= Trt/ dist= binomial link= logit;
RUN;

/**************************************************************
	Section 6.3.2.4
*****************************************************************/
/*fit the gee model with independent patient effect*/
PROC GENMOD data= Dat_L2 descending;
class Trt(ref="1") Time(ref="1") Subject;
model Heal= Trt/ dist= binomial link= logit;
repeated subject= Subject/ corr=ind corrw;
RUN;

/*fit the gee model with independent patient effect*/
PROC GENMOD data= Dat_L2 descending;
class Trt(ref="1") Time(ref="1") Subject;
model Heal= Trt/ dist= binomial link= logit;
repeated subject= Subject/ corr=exch corrw;
RUN;

/**********************************************************
    Section 7.2.2
*******************************************************/
DATA dat6;
input alpha beta; /*alpha is Type-I error and beta is Type-II error*/
pow= 1-beta;
numberator= 2*(quantile("NORMAL",1-alpha/2,0,1)
        +quantile("NORMAL",1-beta,0,1))**2;
datalines;
0.05 0.05
0.05 0.1
0.05 0.15
0.05 0.2
0.05 0.25
0.05 0.3
;
RUN;

/****************************************************************
    Section 7.2.3
    **************************************************/
/*Determine the sample size with two-sided alternative*/
PROC POWER;
twosamplemeans
meandiff= 0.5
stddev= 1
power= 0.8
npergroup= .;
RUN;

/*Determine the sample size with one-sided alternative*/
PROC POWER;
twosamplemeans sides= 1
meandiff= 0.5
stddev= 1
power= 0.8
npergroup= .;
RUN;

/*calculate the statistical power*/
PROC POWER;
twosamplemeans
meandiff= 0.5
stddev= 1
npergroup= 64
power=.;
RUN;

/*calculate the minimum detectable treatment difference*/
PROC POWER;
twosamplemeans
stddev= 1
npergroup= 64
power= 0.8
meandiff= .;
ods output Output= pwout;
RUN;

/*relationship between sample size and statistical power*/
PROC POWER;
twosamplemeans
meandiff= 0.5
stddev= 1
power= 0.2 to 0.9 by 0.05         /*use power from 0.2 to 0.9 by 0.05*/
npergroup= .;
plot x= power min=0.2 max=0.9;
RUN;


/*****************************************************************
    Section 7.2.4
    **********************************************************/
/*equal group size*/
PROC POWER;
twosamplemeans
meandiff= 0.8
stddev= 0.83
power= 0.8
npergroup= .;
RUN;

/*two groups on the order of 2 to 1 ratio*/
PROC POWER;
twosamplemeans
meandiff= 0.8
stddev= 0.83
power= 0.8
groupweights= (2 1)
ntotal= .;
RUN;

/*calculate sample size for unequal variance */
PROC POWER;
twosamplemeans test=diff_satt
meandiff= 2
groupstddevs= 1|2
power= 0.8
ntotal= .
groupweights =(1 2);
RUN;

/*calculate sample size for unequal variance using Wilcoxon test!!!*/
PROC POWER;
twosamplewilcoxon test=wmw
alpha= 0.05
vardist("myordinal1")= ordinal ((0 1 2) : (0.66 0.15 0.19))
vardist("myordinal2")= ordinal ((0 1 2) : (0.61 0.23 0.16))
variables= "myordinal1"|"myordinal2"
power= 0.8
ntotal= .;
RUN;

/******************************************************************
    Section 7.3.1
    **************************************************************/
/*calculate sample size for two independent proportion*/
PROC POWER;
twosamplefreq test=pchi
groupproportions= (0.75 0.5)
power= 0.8
npergroup= .;
RUN;

/*calculate power for two independent proportion*/
PROC POWER;
twosamplefreq test=pchi
groupproportions= (0.75 0.5)
npergroup= 60
power= .;
RUN;

/*relationship between power and significance level*/
PROC POWER;
twosamplefreq test=pchi
groupproportions= (0.75 0.5)
npergroup= 60
alpha= 0.02 to 0.12 by 0.02
power=. ;
RUN;

/*******************************************************************
    Section 7.3.3
    **************************************************************/
/*calculate sample size*/
PROC POWER;
twosamplefreq test=FM
groupproportions= (0.75 0.5)
power= 0.8
npergroup= .;
RUN;

/*calculate sample size with unbalanced group*/
PROC POWER;
twosamplefreq test=FM
sides=1
groupproportions= (0.35 0.2)
groupweights= (2 1)
power= 0.8
ntotal= .;
RUN;

PROC POWER;
twosamplefreq test=FM
sides=1
groupproportions= (0.35 0.2)
power= 0.8
ntotal= .;
RUN;


/***************************************************************
    Section 7.4
    *******************************************************/
PROC POWER;
twosamplesurvival test=logrank
sides=1
gexphs= 0.3|0.2
grouplossexphazards= (0.1 0.1)
accrualtime = 1
followuptime= 2
alpha=0.025
ntotal=.
/*or nevent= .*/
power = 0.8;
RUN;



/*******************************************************************
    Section 7.6.1.4
    ****************************************************************/
Data Long_Trials;
delta= 1.2;         /*minimum treatment effect*/
n= 5;
n1= 2;              /*number of repeated measurement*/
tau= 2;             /*trial duration*/
sig2within= 7;      /*within-subject variability*/
sig2between= 2;     /*between subject variablity*/
alpha= 0.05;        /*significance level*/
pow= 0.9;           /*desired power*/
S1= 40;
sig2= 12*(n-1)*sig2within/(tau**2*n*(n+1))+sig2between;  /*variance of slope*/
sig2_1= 12*(n1-1)*sig2within/(tau**2*n1*(n1+1))+sig2between;
S= ceil((quantile("NORMAL",1-alpha/2,0,1)+quantile("NORMAL",pow,0,1))**2
*2*sig2/delta**2);  /*calculate sample size*/
pow_get= probnorm(sqrt(S1*delta**2/(2*sig2_1))-quantile("NORMAL",1-alpha/2,0,1));
RUN;

/*print result*/
PROC PRINT data= Long_Trials;
var sig2 S pow_get;
RUN;

/****************************************************************************
    Section 7.6.2.2
    *****************************************************************/
Data Long_Trials1;
delta= 1.2;
tau= 2;
sig2within= 7;
sig2between= 2;
alpha= 0.05;
do S= 20 to 100 by 20;
do n= 2 to 10 by 2;
sig2= 12*(n-1)*sig2within/(tau**2*n*(n+1))+sig2between;
pow_get= round(probnorm(sqrt(S*delta**2/(2*sig2))
-quantile("NORMAL",1-alpha/2,0,1)),0.001);
output;
end;
end;
RUN;

/*print the power matrix*/
PROC Transpose data= Long_Trials1;
by S;
id n;
var pow_get;
RUN;

/*******************************************************************
    Section 7.6.2.2
    ****************************************************************/
DATA Long_Trials2;
delta= 0.075;        /*the treatment effect*/
n= 5;               /*number of repeated measurement*/
tau= 2;             /*trial duration*/
rho=0.5;             /*correlation*/
alpha= 0.05;        /*significance level*/
beta=0.1;           /*Type-II error*/
pow= 1-beta;        /*associated power*/
sig2= 3*(n-1)*(1-rho)/(tau**2*n*(n+1));
S= round((quantile("NORMAL",1-alpha/2,0,1)+quantile("NORMAL",pow,0,1))**2
*2*sig2/delta**2);
RUN;

/*print result*/
PROC PRINT data= Long_Trials2;
var sig2 S;
RUN;

/**********************************************************
    Section 7.7.3
    ************************************************/
DATA Long_Trials3;
do CV= 0.1 to 0.7 by 0.1;               /*cv range*/
do PC= 0.1 to 0.5 by 0.05;              /*pv range*/
S= ceil(8*CV**2/PC**2*(1+(1-PC)**2));   /*sample size calculation*/
output;
end;
end;
RUN;

/*print out the result*/
PROC TRANSPOSE data= Long_Trials3;
by CV;
id PC;
var S;
RUN;


/******************************************************
	Section 8.3.3: Amlodipine Data (Continuous data)
*******************************************************/
/* Section 8.3.3.1*/
Data Angina;
infile "Your file path/angina.csv" delimiter="," firstobs=2;
input Protocol nE meanE varE  nC meanC varC;
RUN;
proc iml;
* Read data into IML ;
  use angina;
  read all ;

/* Protocol-specific Analysis*/
/* mean difference*/
iMD = meanE-meanC;
/* var for iMD*/
ivar4MD = (varE/nE+varC/nC);
/* 95% CI for each Study*/
ilowCI = iMD-1.96*sqrt(ivar4MD);
iupCI  = iMD+1.96*sqrt(ivar4MD);
/* inverse variance for weight Calculation */
iwt = 1/ivar4MD;
pctwt = iwt/sum(iwt);
print "Protocol-Specific Summary";
iout = Protocol||iMD||ilowCI||iupCI||pctwt;
print iout;

/*Section 8.3.3.2. Fixed-Effects Meta-Analysis*/
/* MA estimate */
MD = iMD`*pctwt;
/* variance of MD*/
varMD = 1/sum(iwt);
lowCI = MD-1.96*sqrt(varMD);
upCI  = MD+1.96*sqrt(varMD);
/* z-value and p-value*/
z = MD/sqrt(varMD);
pval = 2*(1-probnorm(z));
print "Summary of Fixed-Effects MA";
MA4FE = MD||lowCI||upCI||z||pval;
print MA4FE;

/* Section 8.3.3.3: Random-Effects MA*/
/* calculate the Q, U and Tau-sq*/
C = 8; /* 8 studies*/
U = sum(iwt) - sum(iwt`*iwt)/sum(iwt);
Q = sum(iwt`*(iMD-MD)##2);print Q;
tau2 = (Q-(C-1))/U;
if Q <= C-1 then tau2=0; print tau2;
/*  with tau2, recalculate the weights*/
REiwt = 1/(ivar4MD+tau2);
REpctwt = REiwt/sum(REiwt);print REpctwt;
/* MA estimate */
REMD = iMD`*REpctwt;
/* variance of MD*/
REvarMD = 1/sum(REiwt);
RElowCI = REMD-1.96*sqrt(REvarMD);
REupCI  = REMD+1.96*sqrt(REvarMD);
/* z-value and p-value*/
z = REMD/sqrt(REvarMD);
pval = 2*(1-probnorm(z));
print "Summary of Random-Effects MA";
Re4FE = REMD||RElowCI||REupCI||z||pval;
print RE4FE;


/***************************************************************
Section 9.3.1: Normal-Normal Model
***************************************************************/
/*simulated data*/
DATA S1(keep= y);
n= 30; sigma= 2; mu= 3;
call streaminit(123);
do i= 1 to 30;
y= mu+sigma*rand("Normal");
output;
end;
RUN;

/*the mean and variance of simulated data*/
PROC MEANS data=S1 mean var stddev;
var y;
output out=I1 mean= mean var=var stddev=sd;
RUN;

DATA S2(keep= bayes_norm2norm);
set I1;
n= 30; sigma= 2; mu= 3;
/*the prior parameter*/
mu0= 2; tau0= 0.5;
/*the weight*/
w= tau0**2/(tau0**2+sigma**2/n);
/*the posterior mean*/
muP =w*mean+(1-w)*mu0;
/*the posterior standard deviation*/
sigmaP =sqrt(1/(1/tau0**2+n/sigma**2));
/*set seed*/
call streaminit(123);
/*direct simulation of posterior normal*/
do i= 1 to 10000;
bayes_norm2norm= muP+sigmaP*rand("Normal");
output;
end;
RUN;

/*the quantiles of this posterior distribution*/
PROC UNIVARIATE data= S2;
var bayes_norm2norm;
output out= DS PCTLPTS= 2.5 25 50 75 97.5 PCTLPRE= p;
RUN; 


PROC PRINT data=DS; 
title "Posterior Quantile from Direct Simulation";
RUN;

/*direct simulate using proc*/
PROC MCMC data= S1 outpost= MCsimu seed= 123 nmc= 10000 statistics= all;
parms mu0;
prior mu0 ~ normal(2,sd=0.5);
model y ~ normal(mu0, sd=2);
RUN;

PROC PRINT data=MCsimu;
RUN;

/*the quantiles of this posterior distribution*/
PROC UNIVARIATE data= MCSimu;
var mu0;
output out= MCq PCTLPTS= 2.5 25 50 75 97.5 PCTLPRE = q;
RUN; 

PROC PRINT data=MCq;
RUN;


/***************************************************************
9.3.2 Beta-Binomial Model
***************************************************************/
DATA S4;
input n x;
p= x/n;
datalines;
168 69
182 113
165 120
188 145
RUN;

/*search the root with initial value at (3,3)*/
PROC NLP tech=nmsimp;
min y;
parms a= 3,b= 3;
y=(probbeta(0.5,a,b)-0.75)**2+(probbeta(0.95,a,b)-0.85)**2;
RUN;

/*check the optimization results*/
DATA check;
c1=probbeta(0.5,0.062366,0.183202);
c2=probbeta(0.95,0.062366,0.183202);
RUN;

/*direct simulation*/
DATA S5;
call streaminit(123);
do i= 1 to 10000;
bayes1_betabin= rand("Beta",120.062,45.183);
output;
end;
RUN;

/*print the quantiles*/
PROC UNIVARIATE data= S5;
var bayes1_betabin;
output out= q PCTLPTS= 2.5 25 50 75 97.5 PCTLPRE = q;
RUN; 


/***************************************************************
Section 9.4.1: Blood Pressure Data - Bayesian Linear Regression
***************************************************************/
Data DBP;
infile "Your file path/DBP.csv" delimiter="," firstobs=2; 
input Subject TRT$ DBP1 DBP2 DBP3 DBP4 DBP5 Age Sex$;
diff= DBP5-DBP1;                                         
RUN;

/*fit the Bayes regression model with 1000 burn-in*/
PROC GENMOD data= DBP;
class TRT(ref= "A");
model diff= TRT Age / dist= normal;
bayes nbi=1000 seed= 123 outpost= postDBP ;
RUN;

/***************************************************************
Section 9.4.2: Binomial Data - Bayesian Logistic Regression
***************************************************************/
DATA Betablocker(where= (Center= 1));
infile "Your file path/betablocker.csv" delimiter="," firstobs=2;  
input Deaths Total Center Treatment$;                                       
RUN;

/*fit logistic regression*/
PROC GENMOD data= Betablocker;
class Treatment(ref="Control");
model Deaths/Total= Treatment/dist= binomial link= logit;
RUN;

/*fit the Bayes regression model with 1000 burn-in*/
PROC GENMOD data= BetaBlocker;
class Treatment(ref="Control");
model Deaths/Total= Treatment/dist= binomial link= logit;
bayes nbi=1000 seed= 123 outpost= postbetablocker ;
RUN;


/*fit the Bayes regression model with multivariate normal prior*/
DATA prior;
input _type_ $ Intercept x;
datalines;
Var 1000 1000
Mean 0 0 
;
RUN;

PROC GENMOD data= Betablocker;
class Treatment(ref="Control");
model Deaths/Total= Treatment/dist= binomial link= logit;
bayes nbi=1000 seed= 123 outpost= postbetablocker 
		coeffprior=normal(input= prior);
RUN;


/***************************************************************
Section 9.4.3:Count Data - Bayesian Possion Regression
***************************************************************/
DATA Polyps;
infile "/home/ln23040/polyps.csv" delimiter="," firstobs=2;  
input number treat$ age;                                        
RUN;

/*fit Bayes poisson regression*/
PROC GENMOD data= Polyps;
class treat(ref= "placebo");
model number= treat age / dist= poisson;
bayes nbi=1000 seed= 123 outpost=postPolyps;
RUN;


/*print the quantiles*/
PROC UNIVARIATE data= postPolyps;
var Intercept Treatdrug Age;
output out= q  PCTLPRE=intercept treattdrug age 
	PCTLPTS= 2.5 25 50 75 97.5 PCTLPRE = q;
RUN; 


/****************************************************
	Section 11.3.1. Data manipulation
******************************************************/
Data datAE;
infile "Your file path/datAE.csv" delimiter="," firstobs=2;
input Stage ND2 fD2 ND1 fD1 NC fC;;
proc print; RUN;
proc iml;
* Read data into IML ;
  use datAE;
  read all ;

/* the AE rate*/
p2 = fD2/Nd2; p1 = fD1/ND1; pC=fC/NC;
/* the variance*/
V2 = p2#(1-p2)#(1/ND2);
V1 = p1#(1-p1)#(1/ND1);
VC = pC#(1-pC)#(1/NC);
len=8;

/**********************************************************
	Section 11.3.2. CI method
***********************************************************/
/* just for D1 to Control, same for D2 to Control*/
alpha     = 0.025;
zalpha    = quantile('NORMAL',1-alpha);
lowCI     = (p1-pC)-zalpha*sqrt(v1+vC);
upCI      = (p1-pC)+zalpha*sqrt(v1+vC);
testp     = p1>pC;
diffD1toC = p1-pC;
All       = testp||round(lowCI,0.001)||round(upCI,0.001)
			||round(diffD1toC,0.001)||round(p1,0.001)
			||round(pC,0.001)||ND1||NC;
print "Direct Comparison for D1 to Control";
print All;


/********************************************************
	Section 11.3.3. Indirect CI
*********************************************************/
/* 1. Using normal Approxition */
/*  Calculate CI for Control*/
ClowCI   = pC-zalpha*sqrt(vC);
CupCI    = pC+zalpha*sqrt(vC);
/* CI for Dose 1 */
d1lowCI   = p1-zalpha*sqrt(v1);
d1upCI    = p1+zalpha*sqrt(v1);
/* CI for Dose 2 */
d2lowCI   = p2-zalpha*sqrt(v2);
d2upCI    = p2+zalpha*sqrt(v2);
/* Indirect test for d1/d2 to Control*/
Indirecttest_d1toC = d1lowCI - CupCI;
Indirecttest_d2toC = d2lowCI - CupCI;
print Indirecttest_d1toC;
print Indirecttest_d2toC;

/* 2. Using Exact binomial */
/*  Calculate CI for Control*/
ClowCI   = quantile('BINOM',alpha,pC,nC);
CupCI    = quantile('BINOM',1-alpha,pC,nC);
/* CI for Dose 1 */
d1lowCI   = quantile('BINOM',alpha,p1,nD1);
d1upCI    = quantile('BINOM',1-alpha,p1,nD1);
/* CI for Dose 2 */
d2lowCI   = quantile('BINOM',alpha,p2,nD2);
d2upCI    = quantile('BINOM',1-alpha,p2,nD2);
/* Indirect test for d1/d2 to Control*/
Exact_d1toC = d1lowCI - CupCI>0;
Exact_d2toC = d2lowCI - CupCI>0;
print Exact_d1toC;
print Exact_d2toC;

/* Conclusion for both method: Stage 3*/


/****************************************************
	Section 11.3.4
	Only for "SLM with Normal Approximation"
***************************************************/
proc model data=datAE;
/* the AE rate*/
	p2 = fD2/Nd2; pC=fC/NC;
	/* the variance*/
	V2 = p2*(1-p2)/ND2;
	VC = pC*(1-pC)/NC;
	/* the z-values */
	zalpha = quantile('NORMAL',alpha);
	/* low CI for control */
	ClowCI = (pC-zalpha*sqrt(vC));
	/* upper CI for dose 2 */
	d2upCI = (p2+zalpha*sqrt(v2)); 	
	/* solve the diff =0*/
	eq.diff   = ClowCI-d2upCI;
solve alpha/ out=sol4Norm solveprint;
id fc nc fd2 nd2;
run;

proc print data=sol4Norm; run;







