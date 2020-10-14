/* 연습문제 5.5 */

data melanoma;
input time censor freq @@ ;
cards ;
0 1 9 0 0 0
1 1 6 1 0 1
2 1 2 2 0 4
3 1 1 3 0 5
4 1 2 4 0 3
5 1 0 5 0 17
;
run ;

proc lifetest data=melanoma method=life 
intervals=0 1 2 3 4 5 6 plots=(S,H) graphics ;
time time*censor(0) ;
freq freq ;
title "Life table of 50 skin melanoma patients" ;
run ;


/* 연습문제 5.7 */

proc import datafile="/folders/myfolders/pbc.csv"
out=pbc dbms=csv replace ;
getnames=yes ;
run ;

proc lifetest data=pbc method=KM plots=(S,H) graphics outsurv=pbcsdf ;
time time*censor(0) ;
title "Kaplan-Meier estimation for pbc data" ;
run ; 


/* 연습문제 5.8 */

proc import datafile="/folders/myfolders/gehan.csv"
out=gehan dbms=csv replace ;
getnames=yes ;
run ;

proc lifetest data=gehan method=KM plots=(S,H) graphics outsurv=gehansdf ;
time time*cens(0) ;
title "Kaplan-Meier estimation for gehan data" ;
run ; 

proc lifetest data=gehan method=KM plots=(S,H) graphics outsurv=gehansdf1 ;
time time*cens(0) ;
strata treat ;
title "Kaplan-Meier estimation for gehan data (group : treat)" ;
run ; 

/* 연습문제 5.9 */

data patient ;
input time censor @@ ;
cards ;
11 1 12 1 15 1 44 0 45 1
28 0 16 0 17 0 19 0 30 0
;
run ;

proc lifetest data=patient method=KM plots=(S,H) graphics outsurv=patientsdf ;
time time*censor(0) ;
title "Kaplan-Meier estimation for patient data" ;
run ; 

proc lifetest data=patient method=FH plots=(S,H) graphics outsurv=patientsdf ;
time time*censor(0) ;
title "Nelson-Aalen estimation for patient data" ;
run ;

proc lifetest data=patient method=KM plots=S(ATRISK CL) graphics outsurv=patientsdf ;
time time*censor(0) ;
title "Kaplan-Meier estimation for patient data" ;
run ; 

proc lifetest data=patient method=FH plots=S(ATRISK CL) graphics outsurv=patientsdf ;
time time*censor(0) ;
title "Nelson-Aalen estimation for patient data" ;
run ;












