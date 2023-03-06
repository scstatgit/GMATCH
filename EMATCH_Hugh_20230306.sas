/*Original SAS code: https://support.sas.com/resources/papers/proceedings/proceedings/sugi29/173-29.pdf*/


* Sample data, edited from a SAS example;
* Split them into two datasets for this example.;
data study control;
infile cards;
 rand_num=uniform(3);
 input id study age lwt race smoke ptd ht ui @@;
 if study=1 then output study;
 else output control;
cards;
1 0 14 135 1 0 0 0 0 101 0 14 101 3 1 1 0 0
2 0 15 98 2 0 0 0 0 102 0 15 115 3 0 0 0 1
3 0 16 95 3 0 0 0 0 103 0 16 130 3 0 0 0 0
4 1 17 103 3 0 0 0 0 104 0 17 130 3 1 1 0 1
5 0 17 122 1 1 0 0 0 105 0 17 110 1 1 0 0 0
6 0 17 113 2 0 0 0 0 106 0 17 120 1 1 0 0 0
7 0 17 113 2 0 0 0 0 107 0 17 120 2 0 0 0 0

9 0 18 100 1 1 0 0 0 109 0 18 148 3 0 0 0 0
10 0 18 90 1 1 0 0 1 110 0 18 110 2 1 1 0 0
11 1 19 150 1 0 0 0 0 111 0 19 91 1 1 1 0 1
12 0 19 115 3 0 0 0 0 112 0 19 102 1 0 0 0 0
13 0 19 235 1 1 0 1 0 113 0 19 112 1 1 0 0 1
14 1 20 120 3 0 0 0 1 114 0 20 150 1 1 0 0 0
15 0 20 103 3 0 0 0 0 115 0 20 125 3 0 0 0 1
16 0 20 169 3 0 1 0 1 116 0 20 120 2 1 0 0 0
17 0 20 141 1 0 1 0 1 117 0 20 80 3 1 0 0 1
18 0 20 121 2 1 0 0 0 118 0 20 109 3 0 0 0 0
19 0 20 127 3 0 0 0 0 119 0 20 121 1 1 1 0 1
20 0 20 120 3 0 0 0 0 120 0 20 122 2 1 0 0 0
21 0 20 158 1 0 0 0 0 121 0 20 105 3 0 0 0 0
22 1 21 108 2 1 0 0 1 122 0 21 165 1 1 0 1 0
23 0 21 124 3 0 0 0 0 123 0 21 200 2 0 0 0 0
24 0 21 185 2 1 0 0 0 124 0 21 103 3 0 0 0 0
25 0 21 160 1 0 0 0 0 125 0 21 100 3 0 1 0 0
26 0 21 115 1 0 0 0 0 126 0 21 130 1 1 0 1 0
27 0 22 95 3 0 0 1 0 127 0 22 130 1 1 0 0 0
28 0 22 158 2 0 1 0 0 128 0 22 130 1 1 1 0 1
29 1 23 130 3 0 0 0 0 129 0 23 97 3 0 0 0 1
30 0 23 128 3 0 0 0 0 130 0 23 187 2 1 0 0 0
31 0 23 119 3 0 0 0 0 131 0 23 120 3 0 0 0 0
32 0 23 115 3 1 0 0 0 132 0 23 110 1 1 1 0 0
33 0 23 190 1 0 0 0 0 133 0 23 94 3 1 0 0 0
34 1 24 90 1 1 1 0 0 134 0 24 128 2 0 1 0 0
35 0 24 115 1 0 0 0 0 135 0 24 132 3 0 0 1 0
36 0 24 110 3 0 0 0 0 136 0 24 155 1 1 1 0 0
37 0 24 115 3 0 0 0 0 137 0 24 138 1 0 0 0 0
38 0 24 110 3 0 1 0 0 138 0 24 105 2 1 0 0 0
39 1 25 118 1 1 0 0 0 139 0 25 105 3 0 1 1 0
40 0 25 120 3 0 0 0 1 140 0 25 85 3 0 0 0 1
41 0 25 155 1 0 0 0 0 141 0 25 115 3 0 0 0 0
42 0 25 125 2 0 0 0 0 142 0 25 92 1 1 0 0 0
43 0 25 140 1 0 0 0 0 143 0 25 89 3 0 1 0 0
44 0 25 241 2 0 0 1 0 144 0 25 105 3 0 1 0 0
45 1 26 113 1 1 0 0 0 145 0 26 117 1 1 1 0 0
46 0 26 168 2 1 0 0 0 146 0 26 96 3 0 0 0 0
47 0 26 133 3 1 1 0 0 147 0 26 154 3 0 1 1 0
48 0 26 160 3 0 0 0 0 148 0 26 190 1 1 0 0 0
49 0 27 124 1 1 0 0 0 149 0 27 130 2 0 0 0 1
50 0 28 120 3 0 0 0 0 150 0 28 120 3 1 1 0 1
51 0 28 130 3 0 0 0 0 151 0 28 95 1 1 0 0 0
52 0 29 135 1 0 0 0 0 152 0 29 130 1 0 0 0 1
53 0 30 95 1 1 0 0 0 153 0 30 142 1 1 1 0 0
54 0 31 215 1 1 0 0 0 154 0 31 102 1 1 1 0 0
55 0 32 121 3 0 0 0 0 155 0 32 105 1 1 0 0 0
56 0 34 170 1 0 1 0 0 156 0 34 187 2 1 0 1 0

57 1 26 113 1 1 0 0 0 157 0 26 117 1 1 1 0 0
;
run;
/*8 0 17 119 3 0 0 0 0 108 0 17 142 2 0 0 1 0*/

proc sql;
create table controls_id as 
select
	one.ID as study_id,
	two.ID as control_id,
	one.age as study_age,
	two.age as control_age,
	one.race as study_race,
	two.race as control_race,
	one.rand_num as rand_num
from study one, control two
where (one.age=two.age and one.race=two.race)
;
quit;
* Remove duplicate control subjects;
/*Randomly remove duplicate 3 control subjects*/
proc sort data=controls_id;
 by control_id rand_num;
 run;
 data controls_id;
 set controls_id;
 by control_id rand_num;
 if first.control_id;
 run;
 proc sort data=controls_id ;
 by study_id rand_num;
 run;
data controls_id2 not_enough;
set controls_id;
 by study_id ;
 retain num;
 if first.study_id then num=1;
 if num le 2 then do;
 output controls_id2;
 num=num+1;
 end;
 if last.study_id then do;
 if num le 2 then output not_enough;
 end;
 run;
proc print data=controls_id2(obs=40);
 title2 'matched patients';
 run;
data controls_id3;
merge controls_id2
not_enough(in=b_);
by study_id;
if b_ then delete;
run;
