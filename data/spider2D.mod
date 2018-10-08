param n_blocks := 20;
param M := 10;
param m := 2;
param rhs := M + 2*m + 3;

var x{1..2*n_blocks} >= -1, <= 3;

con01: (M*x[1] +   x[1+1]) + m*x[1+2] +   x[1+3] = rhs - m - 1;
con02: (  x[1] + M*x[1+1]) +   x[1+2] + m*x[1+3] = rhs - m - 1;
con03: m*x[3-2] +   x[3-1] + (M*x[3] +   x[3+1]) + m*x[3+2] +   x[3+3] = rhs;
con04:   x[3-2] + m*x[3-1] + (  x[3] + M*x[3+1]) +   x[3+2] + m*x[3+3] = rhs;
con05: m*x[5-2] +   x[5-1] + (M*x[5] +   x[5+1]) + m*x[5+2] +   x[5+3] = rhs;
con06:   x[5-2] + m*x[5-1] + (  x[5] + M*x[5+1]) +   x[5+2] + m*x[5+3] = rhs;
con07: m*x[7-2] +   x[7-1] + (M*x[7] +   x[7+1]) + m*x[7+2] +   x[7+3] = rhs;
con08:   x[7-2] + m*x[7-1] + (  x[7] + M*x[7+1]) +   x[7+2] + m*x[7+3] = rhs;
con09: m*x[9-2] +   x[9-1] + (M*x[9] +   x[9+1]) + m*x[9+2] +   x[9+3] = rhs;
con10:   x[9-2] + m*x[9-1] + (  x[9] + M*x[9+1]) +   x[9+2] + m*x[9+3] = rhs;
con11: m*x[11-2] +   x[11-1] + (M*x[11] +   x[11+1]) + m*x[11+2] +   x[11+3] = rhs;
con12:   x[11-2] + m*x[11-1] + (  x[11] + M*x[11+1]) +   x[11+2] + m*x[11+3] = rhs;
con13: m*x[13-2] +   x[13-1] + (M*x[13] +   x[13+1]) + m*x[13+2] +   x[13+3] = rhs;
con14:   x[13-2] + m*x[13-1] + (  x[13] + M*x[13+1]) +   x[13+2] + m*x[13+3] = rhs;
con15: m*x[15-2] +   x[15-1] + (M*x[15] +   x[15+1]) + m*x[15+2] +   x[15+3] = rhs;
con16:   x[15-2] + m*x[15-1] + (  x[15] + M*x[15+1]) +   x[15+2] + m*x[15+3] = rhs;
con17: m*x[17-2] +   x[17-1] + (M*x[17] +   x[17+1]) + m*x[17+2] +   x[17+3] = rhs;
con18:   x[17-2] + m*x[17-1] + (  x[17] + M*x[17+1]) +   x[17+2] + m*x[17+3] = rhs;
con19: m*x[19-2] +   x[19-1] + (M*x[19] +   x[19+1]) + m*x[19+2] +   x[19+3] = rhs;
con20:   x[19-2] + m*x[19-1] + (  x[19] + M*x[19+1]) +   x[19+2] + m*x[19+3] = rhs;
con21: m*x[21-2] +   x[21-1] + (M*x[21] +   x[21+1]) + m*x[21+2] +   x[21+3] = rhs;
con22:   x[21-2] + m*x[21-1] + (  x[21] + M*x[21+1]) +   x[21+2] + m*x[21+3] = rhs;
con23: m*x[23-2] +   x[23-1] + (M*x[23] +   x[23+1]) + m*x[23+2] +   x[23+3] = rhs;
con24:   x[23-2] + m*x[23-1] + (  x[23] + M*x[23+1]) +   x[23+2] + m*x[23+3] = rhs;
con25: m*x[25-2] +   x[25-1] + (M*x[25] +   x[25+1]) + m*x[25+2] +   x[25+3] = rhs;
con26:   x[25-2] + m*x[25-1] + (  x[25] + M*x[25+1]) +   x[25+2] + m*x[25+3] = rhs;
con27: m*x[27-2] +   x[27-1] + (M*x[27] +   x[27+1]) + m*x[27+2] +   x[27+3] = rhs;
con28:   x[27-2] + m*x[27-1] + (  x[27] + M*x[27+1]) +   x[27+2] + m*x[27+3] = rhs;
con29: m*x[29-2] +   x[29-1] + (M*x[29] +   x[29+1]) + m*x[29+2] +   x[29+3] = rhs;
con30:   x[29-2] + m*x[29-1] + (  x[29] + M*x[29+1]) +   x[29+2] + m*x[29+3] = rhs;
con31: m*x[31-2] +   x[31-1] + (M*x[31] +   x[31+1]) + m*x[31+2] +   x[31+3] = rhs;
con32:   x[31-2] + m*x[31-1] + (  x[31] + M*x[31+1]) +   x[31+2] + m*x[31+3] = rhs;
con33: m*x[33-2] +   x[33-1] + (M*x[33] +   x[33+1]) + m*x[33+2] +   x[33+3] = rhs;
con34:   x[33-2] + m*x[33-1] + (  x[33] + M*x[33+1]) +   x[33+2] + m*x[33+3] = rhs;
con35: m*x[35-2] +   x[35-1] + (M*x[35] +   x[35+1]) + m*x[35+2] +   x[35+3] = rhs;
con36:   x[35-2] + m*x[35-1] + (  x[35] + M*x[35+1]) +   x[35+2] + m*x[35+3] = rhs;
con37: m*x[37-2] +   x[37-1] + (M*x[37] +   x[37+1]) + m*x[37+2] +   x[37+3] = rhs;
con38:   x[37-2] + m*x[37-1] + (  x[37] + M*x[37+1]) +   x[37+2] + m*x[37+3] = rhs;
con39: m*x[39-2] +   x[39-1] + (M*x[39] +   x[39+1]) = rhs - m - 1;
con40:   x[39-2] + m*x[39-1] + (  x[39] + M*x[39+1]) = rhs - m - 1;

##########################################

option show_stats  1;
option presolve   10;
option substout    1;
option var_bounds  0;
option nl_comments 0;
option nl_permute  0;
option display_precision 0;

option auxfiles rc;

##########################################

print "Assigning suffixes";

suffix blockid IN, integer;

let x[1].blockid := 0+1;
let x[2].blockid := 0+1;
let x[3].blockid := 1+1;
let x[4].blockid := 1+1;
let x[5].blockid := 2+1;
let x[6].blockid := 2+1;
let x[7].blockid := 3+1;
let x[8].blockid := 3+1;
let x[9].blockid := 4+1;
let x[10].blockid := 4+1;
let x[11].blockid := 5+1;
let x[12].blockid := 5+1;
let x[13].blockid := 6+1;
let x[14].blockid := 6+1;
let x[15].blockid := 7+1;
let x[16].blockid := 7+1;
let x[17].blockid := 8+1;
let x[18].blockid := 8+1;
let x[19].blockid := 9+1;
let x[20].blockid := 9+1;
let x[21].blockid := 10+1;
let x[22].blockid := 10+1;
let x[23].blockid := 11+1;
let x[24].blockid := 11+1;
let x[25].blockid := 12+1;
let x[26].blockid := 12+1;
let x[27].blockid := 13+1;
let x[28].blockid := 13+1;
let x[29].blockid := 14+1;
let x[30].blockid := 14+1;
let x[31].blockid := 15+1;
let x[32].blockid := 15+1;
let x[33].blockid := 16+1;
let x[34].blockid := 16+1;
let x[35].blockid := 17+1;
let x[36].blockid := 17+1;
let x[37].blockid := 18+1;
let x[38].blockid := 18+1;
let x[39].blockid := 19+1;
let x[40].blockid := 19+1;

##########################################

let con01.blockid := 0+2;
let con02.blockid := 0+2;
let con03.blockid := 1+2;
let con04.blockid := 1+2;
let con05.blockid := 2+2;
let con06.blockid := 2+2;
let con07.blockid := 3+2;
let con08.blockid := 3+2;
let con09.blockid := 4+2;
let con10.blockid := 4+2;
let con11.blockid := 5+2;
let con12.blockid := 5+2;
let con13.blockid := 6+2;
let con14.blockid := 6+2;
let con15.blockid := 7+2;
let con16.blockid := 7+2;
let con17.blockid := 8+2;
let con18.blockid := 8+2;
let con19.blockid := 9+2;
let con20.blockid := 9+2;
let con21.blockid := 10+2;
let con22.blockid := 10+2;
let con23.blockid := 11+2;
let con24.blockid := 11+2;
let con25.blockid := 12+2;
let con26.blockid := 12+2;
let con27.blockid := 13+2;
let con28.blockid := 13+2;
let con29.blockid := 14+2;
let con30.blockid := 14+2;
let con31.blockid := 15+2;
let con32.blockid := 15+2;
let con33.blockid := 16+2;
let con34.blockid := 16+2;
let con35.blockid := 17+2;
let con36.blockid := 17+2;
let con37.blockid := 18+2;
let con38.blockid := 18+2;
let con39.blockid := 19+2;
let con40.blockid := 19+2;

##########################################

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

write gspider2D;

solve;
display x;
