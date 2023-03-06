/*Original code: https://bioinformaticstools.mayo.edu/research/gmatch/*/


/*------------------------------------------------------------------*
   | The documentation and code below is supplied by HSR CodeXchange.             
   |              
   *------------------------------------------------------------------*/
                                                                                      
                                                                                      
                                                                                      
  /*------------------------------------------------------------------*
   | MACRO NAME  : gmatch
   | SHORT DESC  : Match 1 or more controls to cases using the
   |               GREEDY algorithm
   *------------------------------------------------------------------*
   | CREATED BY  : Kosanke, Jon                  (04/07/2004 16:32)
   |             : Bergstralh, Erik
   *------------------------------------------------------------------*
   | PURPOSE
   |
   | GMATCH Macro to match 1 or more controls for each of N cases
   | using the GREEDY algorithm--REPLACES GREEDY option of MATCH macro.
   | Changes:
   | --cases and controls in same dataset
   | --not mandatory to randomly pre-ort cases and controls, but recommended
   | --options to transform X's and to choose distance metric
   | --input parameters consistent with %DIST macro for optimal matching
   |
   | *******
   |
   | Macro name: %gmatch
   |
   | Authors: Jon Kosanke and Erik Bergstralh
   |
   | Date: July 23, 2003
   |       October 31, 2003...tweaked print/means based on "time" var
   |
   | Macro function:
   |
   | Matching using the GREEDY algorithm
   |
   | The purpose of this macro is to match 1 or more controls(from a total
   | of M) for each of N cases.  The controls may be matched to the cases by
   | one or more factors(X's).  The control selected for a particular
   | case(i) will be the control(j) closest to the case in terms of Dij.
   | Dij can be defined in multiple ways. Common choices are the Euclidean
   | distance and the weighted sum of the absolute differences between the
   | case and control matching factors.  I.e.,
   |
   |     Dij= SQRT [SUM { W.k*(X.ik-X.jk)**2} ],  or
   |
   |     Dij= SUM { W.k*ABS(X.ik-X.jk) },
   |
   |                                      where the sum is over the number
   |                                      of matching factors X(with index
   |                                      k) and W.k = the weight assigned
   |                                      to matching factor k and X.ik =
   |                                      the value of variable X(k) for
   |                                      subject i.
   |
   | The control(j) selected for a case(i) is the one with the smallest Dij
   | (subject to constraints DMAX and DMAXK, defined below). In the case of
   | ties, the first one encountered will be used. The higher the user-defined
   | weight, the more likely it is that the case and control will be matched
   | on the factor.  Assign large weights (relative to the other weights) to
   | obtain exact matches for two-level factors such as gender. An option to
   | using weights might be to standarize the X's in some fashion. The macro
   | has options to standardize all X's to mean 0 and variance 1 and to use
   | ranks.
   |
   | The matching algorithm used is the GREEDY method. Using the greedy method,
   | once a match is made it is never broken.  This may result in inefficiencies
   | if a previously matched control would be a better match for the current
   | case than those controls currently available. (An alternative method is to
   | do optimal matching using the VMATCH & DIST macros. This method guarantees
   | the best possible matched set in terms of minimizing the total Dij.)
   | The GREEDY method generally produces very good matches, especially if the
   | control pool is large relative to the number of cases. When  multiple
   | controls/case are desired, the algorithm first matches 1 control to all
   | cases and then proceeds to select second controls.
   |
   |
   | The gmatch macro checks for missing values of matching variables and the
   | time variable(if specified) and deletes those observations from the input
   | dataset.
   |
   | Call statement:
   |
   |
   | %gmatch(data=,group=,id=,
   |       mvars=,wts=,dmaxk=,dmax=,transf,
   |       time=, dist=,
   |       ncontls=,seedca=,seedco=,
   |       out=,outnmca=,outnmco=,print=);
   |
   | Parameter definitions(R=required parameter):
   |
   |
   |  R    data  SAS data set containing cases and potential controls. Must
   |             contain the ID, GROUP, and the matching variables.
   |
   |  R    group SAS variable defining cases. Group=1 if case, 0 if control.
   |
   |  R     id   SAS CHARACTER ID variable for the cases and controls.
   |
   |
   |  R   mvars  List of numeric matching variables common to both case and
   |             control data sets.  For example, mvars=male age birthyr.
   |
   |  R     wts  List of non-negative weights corresponding to each matching
   |             variable.  For example wts=10 2 1 corresponding to male, age
   |             and birthyr as in the above example.
   |
   |      dmaxk  List of non-negative values corresponding to each matching
   |             variable.  These numbers are the largest possible absolute
   |             differences compatible with a valid match.  Cases will
   |             NOT be matched to a control if ANY of the INDIVIDUAL
   |             matching factor  differences are >DMAXK.  This optional
   |             parameter allows one to form matches of the type male+/-0,
   |             age+/-2, birth year+/-5 by specifying DMAXK=0 2 5.
   |
   |      dmax   Largest value of Dij considered to be a valid match.  If
   |             you want to match exactly on a two-level factor(such as
   |             gender coded as 0 or 1) then assign DMAX to be less than
   |             the weight for the factor.  In the example above, one could
   |             use wt=10 for male and dmax=9.  Leave DMAX blank if any
   |             Dij is a valid match.  One would typically NOT use both
   |             DMAXK and DMAX.  The only advantage to using both, would be
   |             to further restrict potential matches that meet the
   |             DMAXK criteria.
   |
   |       dist  Indicates type of distance to calculate.
   |
   |             1=weighted sum(over matching vars) of
   |             absolute case-control differences(default)
   |
   |             2=weighted Euclidean distance
   |
   |       time  Time variable used for risk set matching.  Matches are only
   |             valid if the control time > case time. May need to
   |
   |     transf  Indicates whether all matching vars are to be transformed
   |             (using the combined case+control data) prior to computing
   |             distances.  0=no(default),
   |                         1=standardize to mean 0 and variance 1,
   |                         2=use ranks of matching variables.
   |
   |    ncontls  Indicates the number of controls to match to each case.  The
   |             default is 1.  With multiple controls per case, the algorithm
   |             will first match every case to one control and then again
   |             match each case to a second control, etc.  Controls selected
   |             on the first pass will be stronger matches than those selected in
   |             later rounds.  The output data set contains a variable (cont_n)
   |             which indicates on which round the control was selected.
   |
   |    seedca   Seed value used to randomly sort the cases prior to
   |             matching. This positive integer will be used as input to
   |             the RANUNI function.  The greedy matching algorithm is
   |             order dependent which, among other things means that
   |             cases matched first will be on average more similar to
   |             their controls than those matched last(as the number of
   |             control choices will be limited).  If the matching order
   |             is related to confounding factors (possibly age or
   |             calendar time) then biases may result.  Therefore it is
   |             generally considered good practice when using the GREEDY
   |             method to randomly sort both the cases and controls
   |             before beginning the matching process.
   |
   |    seedco   Seed value used to randomly sort the controls prior to
   |             matching using the GREEDY method.  This seed value must
   |             also be a positive integer.
   |
   |
   | print= Option to print data for matched cases. Use PRINT=y to
   |        print data and PRINT=n or blank to not print.  Default is y.
   |
   |        out=name of SAS data set containing the results of the matching
   |            process.  Unmatched cases are not included.  See outnm
   |            below.  The default name is __out.  This data set will have
   |            the following layout:
   |
   |          Case_id  Cont_id  Cont_n  Dij  Delta_caco MVARS_ca  MVARS_co
   |             1        67      1     5.2  (Differences & actual
   |             1        78      2     6.1   values for matching factors
   |             2        52      1     2.9   for cases & controls)
   |             2        92      2     3.1
   |             .        .       .      .
   |             .        .       .      .
   |
   |        outnmca=name of SAS data set containing NON-matched cases.
   |                Default name is __nmca .
   |
   |        outnmco=name of SAS data set containing NON-matched controls.
   |                Default name is __nmco .
   |
   |
   |  References:  Bergstralh, EJ and Kosanke JL(1995).  Computerized
   |               matching of controls.  Section of Biostatistics
   |               Technical Report 56.  Mayo Foundation.
   |
   |
   |  Example: 1-1 matching by male(exact), age(+-2) and year(+-5).
   |           The wt for male is not relevant, as only exact matches
   |           on male will be considered.  The weight for age(2) is
   |           double that for year(1).
   |
   |
   |       %gmatch(data=all, group=ca_co,id=clinic,
   |              mvars=male age_od yr_od,
   |              wts=2 2 1, dmaxk=0 2 5,out=mtch,
   |              seedca=87877,seedco=987973);
   |
   *------------------------------------------------------------------*
   | OPERATING SYSTEM COMPATIBILITY
   |
   | UNIX SAS v8   :   YES
   | UNIX SAS v9   :
   | MVS SAS v8    :
   | MVS SAS v9    :
   | PC SAS v8     :
   | PC SAS v9     :
   *------------------------------------------------------------------*
   | EXAMPLES
   |
   | Another example is located at the bottom of the code.
   *------------------------------------------------------------------*
   | Copyright 2004 Mayo Clinic College of Medicine.
   |
   | This program is free software; you can redistribute it and/or
   | modify it under the terms of the GNU General Public License as
   | published by the Free Software Foundation; either version 2 of
   | the License, or (at your option) any later version.
   |
   | This program is distributed in the hope that it will be useful,
   | but WITHOUT ANY WARRANTY; without even the implied warranty of
   | MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   | General Public License for more details.
   *------------------------------------------------------------------*/
 
%MACRO GMATCH(DATA=,GROUP=,ID=,
             MVARS=,WTS=,DMAXK=,DMAX=,DIST=1,
             NCONTLS=1, TIME=,TRANSF=0,
             SEEDCA=,SEEDCO=,PRINT=y,
             OUT=__OUT,OUTNMCA=__NMCA,OUTNMCO=__NMCO);
 
   %LET BAD=0;
 
   %IF %LENGTH(&DATA)=0 %THEN %DO;
      %PUT ERROR: NO DATASET SUPPLIED;
      %LET BAD=1;
   %END;
 
   %IF %LENGTH(&ID)=0 %THEN %DO;
      %PUT ERROR: NO ID VARIABLE SUPPLIED;
      %LET BAD=1;
   %END;
 
   %IF %LENGTH(&GROUP)=0 %THEN %DO;
      %PUT ERROR: NO CASE(1)/CONTROL(0) GROUP VARIABLE SUPPLIED;
      %LET BAD=1;
   %END;
 
   %IF %LENGTH(&MVARS)=0 %THEN %DO;
      %PUT ERROR: NO MATCHING VARIABLES SUPPLIED;
      %LET BAD=1;
   %END;
 
  %IF %LENGTH(&WTS)=0 %THEN %DO;
      %PUT ERROR: NO WEIGHTS SUPPLIED;
      %LET BAD=1;
   %END;
 
   %LET NVAR=0;
   %DO %UNTIL(%SCAN(&MVARS,&NVAR+1,' ')= );
      %LET NVAR=%EVAL(&NVAR+1);
   %END;
   %LET NWTS=0;
   %DO %UNTIL(%QSCAN(&WTS,&NWTS+1,' ')= );
      %LET NWTS=%EVAL(&NWTS+1);
   %END;
   %IF &NVAR^= &NWTS %THEN %DO;
      %PUT ERROR: #VARS MUST EQUAL #WTS;
      %LET BAD=1;
   %END;
 
  %LET NK=0;
   %IF %QUOTE(&DMAXK)^=  %THEN %DO %UNTIL(%QSCAN(&DMAXK,&NK+1,' ')= );
      %LET NK=%EVAL(&NK+1);
   %END;
   %IF &NK>&NVAR %THEN %LET NK=&NVAR;
   %DO I=1 %TO &NVAR;
      %LET V&I=%SCAN(&MVARS,&I,' ');
   %END;
 
  %IF &NWTS>0 %THEN %DO;
        DATA _NULL_;
        %DO I=1 %TO &NWTS;
             %LET W&I=%SCAN(&WTS,&I,' ');
             IF &&W&I<0 THEN DO;
                  PUT 'ERROR: WEIGHTS MUST BE NON-NEGATIVE';
                  CALL SYMPUT('BAD','1');
             END;
        %END;
        RUN;
   %END;
 
  %IF &NK>0 %THEN %DO;
        DATA _NULL_;
        %DO I=1 %TO &NK;
             %LET K&I=%SCAN(&DMAXK,&I,' ');
             IF &&K&I<0 THEN DO;
                  PUT 'ERROR: DMAXK VALUES MUST BE NON-NEGATIVE';
                  CALL SYMPUT('BAD','1');
             END;
        %END;
        RUN;
   %END;
 
    %MACRO MAX1;
      %IF &DMAX^= %THEN %DO;
         & __D<=&DMAX
      %END;
      %DO I=1 %TO &NK;
         & ABS(__CA&I-__CO&I)<=&&K&I
      %END;
    %MEND MAX1;
 
   %macro greedy;
    %GLOBAL BAD2;
 
      data __CHECK; set &DATA;
          __id=&id;
          if __id="" then delete;
          %DO I=1 %TO &NVAR;
                IF %scan(&mvars,&i)=. THEN DELETE;
           %END;
           %IF &TIME^= %THEN %DO;
                IF &TIME=. THEN DELETE;
           %END;
       run;
 
      *** transform data if requested/separate cases & controls;
      %if &transf=1 %then %do;
      proc standard data=__check m=0 s=1 out=_stdzd; var &mvars;
      data _caco;
        set _stdzd;
      %end;
 
      %if &transf=2 %then %do;
      proc rank data=__check out=_ranks; var &mvars;
      data _caco;
        set _ranks;
      %end;
 
      %if &transf=0 %then %do;
      data _caco;
        set __check;
      %end;
 
 
      DATA __CASE; SET _caco;
           if &group=1;
      DATA __CASE; SET __CASE END=EOF;
       KEEP __IDCA __CA1-__CA&NVAR __R &mvars
         %if &time^= %then %do;
             __catime
         %end;
          ;
         __IDCA=&ID;
         %if &time^= %then %do;
            __catime=&time;
         %end;
         %DO I=1 %TO &NVAR;
            __CA&I=&&V&I;
         %END;
         %if &seedca^= %then %do;
         SEED=&SEEDCA;
         __R=RANUNI( SEED  );
         %end;
         %else %do;
         __R=1;
         %end;
 
         IF EOF THEN CALL SYMPUT('NCA',_N_);
      PROC SORT; BY __R __IDCA;
 
      DATA __CONT; SET _caco;
         if &group=0;
      DATA __CONT; SET __CONT END=EOF;
       KEEP __IDCO __CO1-__CO&NVAR __R &mvars
        %if &time^= %then %do;
           __cotime
        %end;
        ;
         __IDCO=&ID;
         %if &time^= %then %do;
            __cotime=&time;
         %end;
         %DO I=1 %TO &NVAR;
            __CO&I=&&V&I;
         %END;
         %if &seedco^= %then %do;
         SEED=&SEEDCo;
         __R=RANUNI( SEED  );
         %end;
         %else %do;
         __R=1;
         %end;
 
         IF EOF THEN CALL SYMPUT('NCO',_N_);
      RUN;
      %LET BAD2=0;
      %IF &NCO < %EVAL(&NCA*&NCONTLS) %THEN %DO;
         %PUT ERROR: NOT ENOUGH CONTROLS TO MAKE REQUESTED MATCHES;
         %LET BAD2=1;
      %END;
 
      %IF &BAD2=0 %THEN %DO;
         PROC SORT; BY __R __IDCO;
         DATA __MATCH;
          KEEP __IDCA __CA1-__CA&NVAR __DIJ __MATCH __CONT_N
          %if &time^= %then %do;
             __catime __cotime
          %end;
          ;
          ARRAY __USED(&NCO) $ 1 _TEMPORARY_;
            DO __I=1 TO &NCO;
               __USED(__I)='0';
            END;
            DO __I=1 TO &NCONTLS;
               DO __J=1 TO &NCA;
                  SET __CASE POINT=__J;
                  __SMALL=.;
                  __MATCH=.;
                  DO __K=1 TO &NCO;
                     IF __USED(__K)='0' THEN DO;
                        SET __CONT POINT=__K;
 
                       %if &dist=2 %then %do;
                        **wtd euclidian dist;
                         __D= sqrt(
                         %do k=1 %to &nvar;
                         %scan(&wts,&k)*(__ca&k - __co&k)**2
                         %if &k<&nvar %then + ;
                        %end;
                         );
                       %end;
                       %else %do;
                        **wtd sum absolute diff;
                         __D=
                        %do k=1 %to &nvar;
                        %scan(&wts,&k)*abs(__ca&k - __co&k )
                        %if &k<&nvar %then + ;
                        %end;
                          ;
                       %end;
 
                        IF __d^=. & (__SMALL=. | __D<__SMALL) %MAX1
                        %if &time^= %then %do;
                           & __cotime > __catime
                        %end;
                        THEN DO;
                           __SMALL=__D;
                           __MATCH=__K;
                           __DIJ=__D;
                           __CONT_N=__I;
                        END;
                     END;
                  END;
                  IF __MATCH^=. THEN DO;
                     __USED(__MATCH)='1';
                     OUTPUT;
                  END;
               END;
            END;
            STOP;
         DATA &OUT;
          SET __MATCH;
          SET __CONT POINT=__MATCH;
          KEEP __IDCA __IDCO __CONT_N __DIJ __CA1-__CA&NVAR
               __CO1-__CO&NVAR __d1-__d&nvar __absd1-__absd&nvar  __WT1-__WT&NVAR
                  __catime __cotime __dtime;
 
          %if &time= %then %do;
              __cotime=.; __catime=.;
          %end;
          LABEL
                   __catime="&time/CASE"
                   __cotime="&time/CONTROL"
                   __dtime="&time/ABS. DIFF"
                __CONT_N='CONTROL/NUMBER'
                __DIJ='DISTANCE/D_IJ'
               %DO I=1 %TO &NVAR;
                __CA&I="&&V&I/CASE"
                __CO&I="&&V&I/CONTROL"
                __absd&I="&&V&I/ABS. DIFF "
                __d&I="&&V&I/DIFF "
                __WT&I="&&V&I/WEIGHT"
              %END;
                ;
             %DO I=1 %TO &NVAR;
                __d&i= (__CA&I-__CO&I);      **raw diff;
                __absd&I=abs(__CA&I-__CO&I); **abs diff;
                __WT&I=&&W&I;
             %END;
                __dtime=__cotime-__catime;
 
         PROC SORT DATA=&OUT; BY __IDCA __CONT_N;
         proc sort data=__case; by __IDCA;
         data &outnmca; merge __case
              &out(in=__inout where=(__cont_n=1)); by __idca;
              if __inout=0; **non-matches;
 
         proc sort data=__cont; by __IDCO;
         proc sort data=&out; by __IDCO;
         data &outnmco; merge __cont
              &out(in=__inout); by __idco;
              if __inout=0; **non-matched controls;
         proc sort data=&out; by __IDCA; **re-sort by case id;
 
       %if %upcase(&print)=Y %then %do;
         PROC PRINT data=&out LABEL SPLIT='/';
          VAR __IDCA __IDCO __CONT_N
 
           __DIJ
          %DO I=1 %TO &NVAR;
           __absd&I
          %END;
          %if &time^= %then %do;
           __dtime
          %end;
          %DO I=1 %TO &NVAR;
           __CA&I __CO&I
          %END;
          %if &time^= %then %do;
           __catime __cotime
          %end;
           ;
          sum __dij;
 
         title9'Data listing for matched cases and controls';
         footnote"Greedy matching(gmatch) macro: data=&data group=&group id=&id    ";
         footnote2"   mvars=&mvars  wts=&wts dmaxk=&dmaxk dmax=&dmax ncontls=&ncontls";
         footnote3"   transf=&transf dist=&dist time=&time seedca=&seedca  seedco=&seedco";
         footnote4"   out=&out   outnmca=&outnmca  outnmco=&outnmco";
         run;
         title9'Summary data for matched cases and controls--one obs/control';
          %if &sysver ge 8 %then %do;
         proc means data=&out  maxdec=3 fw=8
           n mean median min p10 p25 p75 p90 max sum;
         %end;
         %else %do;
         proc means data=&out maxdec=3
          n mean min max sum;
         %end;
         class __cont_n;
          var __dij
 
              %do I=1 %TO &NVAR;
                  __absd&I
              %end;
              %if &time^= %then %do;
                  __dtime
              %end;
              %do I=1 %TO &NVAR;
                  __ca&I
              %end;
              %if &time^= %then %do;
                  __catime
              %end;
              %do I=1 %TO &NVAR;
                  __co&I
              %end;
              %if &time^= %then %do;
                  __cotime
              %end;
                 ;
         run;
         *** estimate matching var means within matched sets for controls;
         proc means data=&out  n mean noprint; by __idca;
          var __dij
         %do i=1 %to &nvar;
            __co&i
         %end;
              __cotime
            ;
         output out=_mcont n=n_co mean=__dijm
         %do i=1 %to &nvar;
           __com&i
         %end;
             __tcom
           ;
         data _onecase; set &out; by __idca; if first.__idca;
         data __camcon; merge _onecase _mcont; by __idca;
 
         keep __idca n_co __dijm
             __dtime __catime  __tcom
          %do i=1 %to &nvar;
           __ca&i __com&i  __actd&i __absd&i
          %end;
         ;
 
 
         %do i=1 %to &nvar;
         __absd&i=abs(__ca&i - __com&i);
         __actd&i=(__ca&i - __com&i);
        %end;
         __dtime=__tcom-__catime
          ;
 
       label
        n_co="No./CONTROLS"
        __dijm="Average/Dij"
        __dtime="&time/Mean Time DIFF"
        __tcom="&time/Mean CONT TIME"
 
       %do i=1 %to &nvar; %let vvar=%scan(&mvars,&i);
         __absd&i="&vvar/Mean ABS. DIFF"
         __com&i="&vvar/Mean CONTROL"
       %end;
         ;
      title9'Summary data for matched cases and controls--one obs/case(using average control value)';
      %if &sysver ge 8 %then %do;
      proc means data=__camcon maxdec=3 fw=8
        n mean median min p10 p25 p75 p90 max sum;
      %end;
      %else %do;
      proc means data=__camcon maxdec=3
        n mean min max sum;
      %end;
      var n_co __dijm
      %do i=1 %to &nvar;
       __absd&i
      %end;
      %if &time^= %then %do;
       __dtime
      %end;
      %do i=1 %to &nvar;
      __ca&i
      %end;
      %if &time^= %then %do;
       __catime
      %end;
      %do i=1 %to &nvar;
      __com&i
      %end;
      %if &time^= %then %do;
      __tcom
      %end;
          ;
    %end; **end of print=y loop**;
   %END;  **end of bad2=0 loop**;
   run;
   title9; footnote;
   run;
 
   %mend greedy;
 
   %IF &BAD=0 %THEN %DO;
         %GREEDY
   %END;
%MEND GMATCH;
 
 
*********************************************************************************;


* Sample data, edited from a SAS example;
* Split them into two datasets for this example.;
data study control;
infile cards;
 rand_num=uniform(0);
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

data tmp;
set study control;
run;


*********************************************************************************;


/*ex) 1:5 matching*/
%let nctr=5;
%gmatch(data=tmp, group=study, id=id,
       mvars=age race,wts=1 1,dmaxk= 0 0, transf=0,
       dist=1, ncontls=&nctr, seedca=1234, seedco=1234, print=Y);
run;


/*remove of 1:(equal or less than N-1) matched case*/
/*in this case, N=5*/
proc sort data=__out; by __idca __cont_n; run;
data __out2 __nenough;
set __out;
by __idca;
retain __cont_n;
if first.__idca then __cont_n=1;
if __cont_n le &nctr then do;
output __out2;
__cont_n=__cont_n+1;
end;
if last.__idca then do;
if __cont_n le &nctr then output __nenough;
end;
run;
data __out3;
merge __out2 __nenough(in=b_);
by __idca;
if b_ then delete;
run;
