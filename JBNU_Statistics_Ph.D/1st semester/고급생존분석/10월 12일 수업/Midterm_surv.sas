data pro5 ;
input time censor @@ ;
cards ;
3 1 4 1 5 0 6 1 6 0
8 0 11 1 14 1 15 1 16 0
;
run ;
proc lifereg data=pro5 ;
model time*censor(0)= / dist=weibull ;
run ;
proc lifetest data=pro5 method=KM plots=S(ATRISK CL) graphics ;
time time*censor(0) ;
title "Kaplan-Meier estimation for patient data" ;
run ; 
proc lifetest data=pro5 method=FH plots=S(ATRISK CL) graphics ;
time time*censor(0) ;
title "Nelson-Aalen estimation for patient data" ;
run ;