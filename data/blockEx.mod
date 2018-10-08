param n:=20;
param M:=10;

param rhs{1..n} default M+2;
let rhs[1]:=M+1;
let rhs[n]:=M+1;

var x{0..n+1} >= -1, <= 3;

def_beg: 
  x[0]   = 0;
def_end: 
  x[n+1] = 0;

con{i in 1..n}: 
  x[i-1] + M*x[i] + x[i+1] = rhs[i];
  
##########################################

option show_stats  1;
option presolve   10;
option substout    1;
option var_bounds  0;
option nl_comments 0;
option nl_permute  0;
option display_precision 0;

option auxfiles rc;

################################################################################

print "Assigning suffixes";

suffix blockid IN, integer;
 
for {i in 1..n}
  let x[i].blockid := i;

################################################################################

for {i in 1..n}
  let con[i].blockid := i+1;

################################################################################

print "Before write command";

write gblockEx;

print "Done with write";

########################################################################

print "Constraint suffixes missing: ";

for {k in 1.._sncons} { 
  
  if _scon[k].blockid==0 then
  
    print "  ",_sconname[k];
}

print "Variable suffixes missing: ";

for {k in 1.._snvars} {

  if _svar[k].blockid==0 then
  
    print "  ",_svarname[k];
}

end;

print "===========================================";

print "Constraints: ";

for {k in 1.._sncons}
  print _sconname[k],"  ",_scon[k].blockid;
  
print "Variables: ";

for {k in 1.._snvars}
    print _svarname[k],"  ",_svar[k].blockid;
    
option solver "/home/ali/ampl/ipopt";

solve;

display x;