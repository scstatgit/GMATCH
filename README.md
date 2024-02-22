# MATCH with SAS 9.4


## EMATCH
Exact matching can be used in case-control matching for epidemiological studies. However, I am having trouble with the 'Remove duplicate control group' algorithm in the original SAS code, which can be found at the URL below for exact matching. If you are experiencing the same issue, you can use this code to solve it. </br>
Original SAS code: https://support.sas.com/resources/papers/proceedings/proceedings/sugi29/173-29.pdf <br>
*maybe there is some error in this code. so don't use my uploaded EMATCH.sas

## GMATCH
Greedy match can be used in case-control 1:N matching for epidemiological studies. However, I am having trouble with the 'remove of 1:(equal or less than N-1) matched case' in the original SAS code, which can be found at the URL below for exact matching. If you are experiencing the same issue, you can use this code to solve it. </br>


## PSMatch_Multi
PS match can be used in case-control 1:N matching for epidemiological studies. </br>
This code is PS match with simple example. </br>
Merging matched data on pscore by id with sequence and group. <br>
Original SAS code: https://support.sas.com/resources/papers/proceedings17/0812-2017.pdf


