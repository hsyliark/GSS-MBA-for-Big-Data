/* problem 6.8 */
data ex ;
input x1 x2 treat @@ ;
cards ;
6 7 1 5 9 1 8 6 1 4 9 1
7 9 1 3 3 2 1 6 2 2 3 2
2 3 3 5 1 3 3 1 3 2 3 3
;
run ;
proc glm data=ex ;
class treat ;
model x1 x2 = treat / ss3 ;
manova h = treat / printe ;
means treat ;
run ;

/* problem 6.19 */
proc import datafile="/folders/myfolders/problem 6_19.csv"
out=fuel dbms=csv replace ;
getnames=yes ;
run ;
proc glm data=fuel ;
class type ;
model fuel repair capital = type / ss3 ;
manova h = type / printe ;
means type ;
run ;
proc import datafile="/folders/myfolders/problem 6_19_correct.csv"
out=fuel1 dbms=csv replace ;
getnames=yes ;
run ;
proc glm data=fuel1 ;
class type ;
model fuel repair capital = type / ss3 ;
manova h = type / printe ;
means type ;
run ;

