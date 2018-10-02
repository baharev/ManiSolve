
# NUMBER OF STAGES
param N := 54;

# FEED STAGE LOCATION
param N_F := 30;

# NUMBER OF COMPONENTS
param C := 3;

# DISTILLATE MOLAR FLOW RATE
param D;

# Bad: 459, 465, 455, 45, 44, 435
# Good: 0.434
let D := 0.4342;

# VAPOR FLOW RATE, FROM SPECIFICATION
param V := 1.38;

# FEED LOCATION AND VALUES ARE GIVEN IN THE DATA SECTION

param F{j in 0..N};
param f{i in 1..C, j in 0..N};

# AUXILIARY PARAMETERS

param L{j in 0..N};
param B := F[N_F] - D;

################## PARAMETERS #########################################

param a{1..C};
param b{1..C};
param c{1..C};

param r{1..C, 1..C};
param s{1..C, 1..C};

param P := 100000.0;

# LOWER/UPPER BOUNDS ON THE VARIABLES, INITIAL ESTIMATES (IF NEEDED)
# VALUES ARE GIVEN IN THE DATA SECTION

param x_L{1..C, 1..N}; param x_U{1..C, 1..N}; param x_0{1..C, 1..N};

param T_0{1..N};

param shift = 336.3;

param scale = 383.4 - 336.3;

############### VARIABLES #############################################

var x{i in 1..C, j in 1..N} >= x_L[i,j], <= x_U[i,j], := x_0[i,j];

var T{j in 1..N} >= 0.0, <= 1.0, := T_0[j];

####### DEFINED VARIABLES (THEY ARE ELIMINATED BY PRESOLVE / SUBTITUTION) ######

var Tint{j in 1..N} = scale*T[j] + shift;

var p{i in 1..C, j in 1..N} = exp(a[i]+b[i]/(Tint[j]+c[i]));

var rcp_T{j in 1..N} = 1.0/Tint[j];

var Lambda{i1 in 1..C, i2 in 1..C, j in 1..N} = exp(r[i1,i2]+s[i1,i2]*rcp_T[j]);

var sum_xLambda{i in 1..C, j in 1..N} = sum{i1 in 1..C} (x[i1,j]*Lambda[i,i1,j]);

var rcp_sum_xLambda{i in 1..C, j in 1..N} = 1.0/sum_xLambda[i,j];

var gamma{i in 1..C, j in 1..N} =
  exp( -log(sum_xLambda[i,j]) + 1.0 - (sum{i2 in 1..C} (x[i2,j]*Lambda[i2,i,j]*rcp_sum_xLambda[i2,j])) );

var K{i in 1..C, j in 1..N} = gamma[i,j]*(p[i,j]/P);

var y{i in 1..C, j in 1..N} = K[i,j]*x[i,j];
  
############## EQUATIONS ##############################################

M_eq_N{i in {1, 3}}:
	L[N-1]*x[i,N-1] - B*x[i,N] - V*y[i,N] = 0.0;

M_eq{j in 2..N-1, i in {1, 3}}:
	L[j]*x[i,j] + V*y[i,j] - f[i,j] - L[j-1]*x[i,j-1] - V*y[i,j+1] = 0.0;
	
Eq{j in 1..N}:
    sum{i in 1..C} y[i,j] - 1.0 = 0.0;

S_x_eq{j in 1..N}:
	sum{i in 1..C} x[i,j] - 1.0 = 0.0;

M_eq_1{i in {1, 3}}:
	L[1]*x[i,1] + D*y[i,1] - V*y[i,2] = 0.0;

################### DATA SECTION ######################################

data;

let a[1] := 23.4832;
let a[2] := 20.5110;
let a[3] := 20.9064;

let b[1] := -3634.01;
let b[2] := -2664.30;
let b[3] := -3096.52;

let c[1] := -33.768;
let c[2] := -79.483;
let c[3] := -53.668;

let r[1,2] :=  0.7411;
let r[1,3] :=  0.9645;
let r[2,3] := -1.4350;

let r[2,1] := -1.0250;
let r[3,1] := -0.9645;
let r[3,2] :=  2.7470;

let r[1,1] := 0.0;
let r[2,2] := 0.0;
let r[3,3] := 0.0;

let s[1,2] := -477.00;
let s[1,3] := -903.1024;
let s[2,3] :=  768.20;

let s[2,1] :=  72.78;
let s[3,1] := -140.9995;
let s[3,2] := -1419.0;

let s[1,1] := 0.0;
let s[2,2] := 0.0;
let s[3,3] := 0.0;

# LOWER AND UPPER BOUNDS ON THE VARIABLES

for {j in 1..N} {
	let x_L[1,j] := 0.0;
	let x_U[1,j] := 1.0;
	let x_L[2,j] := 0.0;
	let x_U[2,j] := 1.0;
	let x_L[3,j] := 0.0;
	let x_U[3,j] := 1.0;
}

# THIS BOUND SEEMS TO BE REASONABLE FROM ENGINEENERING POINT OF VIEW

# let x_L[1,1] := 0.83;

# FEED VALUES, LOCATION

for {i in 1..C, j in 0..N}
	let f[i,j] := 0.0;

let f[1, N_F] := 0.4098370;
let f[2, N_F] := 0.01229769;
let f[3, N_F] := 0.06090665;

for {j in 0..N}
	let F[j] := sum{i in 1..C} f[i,j];

for {j in 0..N}
  let L[j] := if (j!=N) then V-D+sum{k in 0..j}F[k] else  -D+sum{k in 0..j}F[k];

#######################################################################

# DUMB INITIAL ESTIMATES (IF NEEDED)

for {i in 1..C, j in 1..N}
	let x_0[i,j] := 0.33;

for {j in 1..N}
	let T_0[j] := 0.001;

#######################################################################

option show_stats 1;
option presolve 10;
option substout 1;
option var_bounds 2;
option nl_comments 1;

#######################################################################

suffix blockid IN, integer;

for {i in 1..C, j in 1..N} {
  let x[i,j].blockid := j+1;
  let T[j].blockid   := j+1;
}

let x[1,1].blockid := 1;
let x[3,1].blockid := 1;

#######################################################################

for {i in {1, 3}}
   let M_eq_1[i].blockid := 3;
   
for {i in  {1, C}, j in 2..N-1}
  let M_eq[j,i].blockid := j+2;

for {i in {1, C}}
   let M_eq_N[i].blockid := N+2;

for {j in 1..N}
  let Eq[j].blockid := j+1;
   
for {i in 1..C, j in 1..N}
  let S_x_eq[j].blockid := j+1;

option auxfiles rc;

write gmss54_B;

#######################################################################

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

#######################################################################

print "Variable blocks: ";

for {blk in 0..N+3} {
    print;
    print "Block", blk;
    for {k in 1.._snvars} {
        if _svar[k].blockid==blk then
            print "  ",_svarname[k];
    }
}

#######################################################################

print "Constraint blocks: ";

for {blk in 0..N+3} {
    print;
    print "Block", blk;
    for {k in 1.._sncons} {
        if _scon[k].blockid==blk then
            print "  ",_sconname[k];
    }
}

#######################################################################

#end;

option solver ipopt;

option seed 31;

param N_TRIALS := 1000;

printf "@@@ N_TRIALS\n%d\n", N_TRIALS;

param N_SOL := N_TRIALS;

param sol{1..N_SOL} default 0;
param val{1..N_SOL} default 0;
param failed default 0;
param selected_component;

solexpand _scon[1];

print "~~~~~";
print _snvars;
for {i in 1.._snvars}
  printf "%s\n", _svarname[i];

print "";
  
for {1..N_TRIALS} {

  for {j in 1.._snvars}
    let _svar[j] := Uniform(_svar[j].lb, _svar[j].ub);

  solve;

  display solve_result, solve_result_num, solve_message, solve_exitcode;

  if (solve_result_num < 200) then {
  
    let selected_component := x[2,ceil(N/3)];
   
    for {i in 1..N_SOL} {
  
      if (val[i]==0.0 or (abs(selected_component-val[i]) < 1.0e-4)) then { 
    
        let val[i] := selected_component;
      
        let sol[i] := sol[i] + 1;
        
        display selected_component;
      
        printf "@@@ %d\n", i;

        for {k in 1.._snvars}
          printf "%.15g\n", _svar[k];        

#  To dump the unscaled temperatures
#        for {k in 1.._snvars-N}
#          printf "%.15g\n", _svar[k];
#          
#        for {k in 1..N}
#          printf "%.15g\n", Tint[k];
      
        print "";
      
        break;
      } # if (val[i]==0.0 ...
      
    } # for {i in 1..N_SOL} 
    
  } # if (solve_result_num < 200) then
  else {
    print "FAILED";
    let failed := failed + 1;
  }
} # for {1..N_TRIALS}

display failed;

printf "frequency of solutions:\n\n";

for {k in 1..N_TRIALS} {

  if (sol[k]!=0) then printf "%d  %d\n", k, sol[k];
}

printf "\n==============================\n";
printf "solutions: %d out of %d trials\n\n", sum {k in 1..N_SOL} sol[k], N_TRIALS;

end;
