OPENQASM 2.0;
include "qelib1.inc";
qreg q[30];
creg c[30];
h q[4];
cx q[3],q[4];
tdg q[4];
cx q[0],q[4];
t q[4];
cx q[3],q[4];
cx q[0],q[4];
cx q[0],q[3];
tdg q[3];
cx q[0],q[3];
s q[4];
h q[5];
cx q[2],q[5];
t q[5];
cx q[1],q[5];
t q[5];
cx q[2],q[5];
cx q[1],q[5];
cx q[1],q[2];
tdg q[2];
cx q[1],q[2];
sdg q[5];
h q[10];
cx q[9],q[10];
tdg q[10];
cx q[6],q[10];
t q[10];
cx q[9],q[10];
cx q[6],q[10];
cx q[6],q[9];
tdg q[9];
cx q[6],q[9];
s q[10];
h q[11];
cx q[8],q[11];
t q[11];
cx q[7],q[11];
t q[11];
cx q[8],q[11];
cx q[7],q[11];
cx q[7],q[8];
tdg q[8];
cx q[7],q[8];
sdg q[11];
h q[18];
cx q[17],q[18];
tdg q[18];
cx q[14],q[18];
t q[18];
cx q[17],q[18];
cx q[14],q[18];
cx q[14],q[17];
tdg q[17];
cx q[14],q[17];
s q[18];
h q[19];
cx q[16],q[19];
t q[19];
cx q[15],q[19];
t q[19];
cx q[16],q[19];
cx q[15],q[19];
cx q[15],q[16];
tdg q[16];
cx q[15],q[16];
sdg q[19];
h q[24];
cx q[23],q[24];
tdg q[24];
cx q[20],q[24];
t q[24];
cx q[23],q[24];
cx q[20],q[24];
cx q[20],q[23];
tdg q[23];
cx q[20],q[23];
s q[24];
h q[25];
cx q[22],q[25];
t q[25];
cx q[21],q[25];
t q[25];
cx q[22],q[25];
cx q[21],q[25];
cx q[21],q[22];
tdg q[22];
cx q[21],q[22];
sdg q[25];
cx q[2],q[4];
tdg q[4];
cx q[0],q[4];
tdg q[4];
cx q[2],q[4];
cx q[0],q[4];
cx q[0],q[2];
t q[2];
cx q[0],q[2];
h q[4];
cx q[3],q[5];
t q[5];
cx q[1],q[5];
tdg q[5];
cx q[3],q[5];
cx q[1],q[5];
cx q[1],q[3];
t q[3];
cx q[1],q[3];
h q[5];
cx q[8],q[10];
tdg q[10];
cx q[6],q[10];
tdg q[10];
cx q[8],q[10];
cx q[6],q[10];
cx q[6],q[8];
t q[8];
cx q[6],q[8];
h q[10];
cx q[9],q[11];
t q[11];
cx q[7],q[11];
tdg q[11];
cx q[9],q[11];
cx q[7],q[11];
cx q[7],q[9];
t q[9];
cx q[7],q[9];
h q[11];
cx q[16],q[18];
tdg q[18];
cx q[14],q[18];
tdg q[18];
cx q[16],q[18];
cx q[14],q[18];
cx q[14],q[16];
t q[16];
cx q[14],q[16];
h q[18];
cx q[17],q[19];
t q[19];
cx q[15],q[19];
tdg q[19];
cx q[17],q[19];
cx q[15],q[19];
cx q[15],q[17];
t q[17];
cx q[15],q[17];
h q[19];
cx q[22],q[24];
tdg q[24];
cx q[20],q[24];
tdg q[24];
cx q[22],q[24];
cx q[20],q[24];
cx q[20],q[22];
t q[22];
cx q[20],q[22];
h q[24];
cx q[23],q[25];
t q[25];
cx q[21],q[25];
tdg q[25];
cx q[23],q[25];
cx q[21],q[25];
cx q[21],q[23];
t q[23];
cx q[21],q[23];
h q[25];
h q[12];
cx q[10],q[12];
t q[12];
cx q[4],q[12];
t q[12];
cx q[10],q[12];
cx q[4],q[12];
cx q[4],q[10];
tdg q[10];
cx q[4],q[10];
sdg q[12];
h q[13];
cx q[11],q[13];
tdg q[13];
cx q[5],q[13];
t q[13];
cx q[11],q[13];
cx q[5],q[13];
cx q[5],q[11];
tdg q[11];
cx q[5],q[11];
s q[13];
h q[26];
cx q[24],q[26];
t q[26];
cx q[18],q[26];
t q[26];
cx q[24],q[26];
cx q[18],q[26];
cx q[18],q[24];
tdg q[24];
cx q[18],q[24];
sdg q[26];
h q[27];
cx q[25],q[27];
tdg q[27];
cx q[19],q[27];
t q[27];
cx q[25],q[27];
cx q[19],q[27];
cx q[19],q[25];
tdg q[25];
cx q[19],q[25];
s q[27];
cx q[10],q[13];
tdg q[13];
cx q[5],q[13];
tdg q[13];
cx q[10],q[13];
cx q[5],q[13];
cx q[5],q[10];
t q[10];
cx q[5],q[10];
h q[13];
cx q[11],q[12];
t q[12];
cx q[4],q[12];
tdg q[12];
cx q[11],q[12];
cx q[4],q[12];
cx q[4],q[11];
t q[11];
cx q[4],q[11];
h q[12];
cx q[24],q[27];
tdg q[27];
cx q[19],q[27];
tdg q[27];
cx q[24],q[27];
cx q[19],q[27];
cx q[19],q[24];
t q[24];
cx q[19],q[24];
h q[27];
cx q[25],q[26];
t q[26];
cx q[18],q[26];
tdg q[26];
cx q[25],q[26];
cx q[18],q[26];
cx q[18],q[25];
t q[25];
cx q[18],q[25];
h q[26];
h q[28];
cx q[26],q[28];
t q[28];
cx q[12],q[28];
t q[28];
cx q[26],q[28];
cx q[12],q[28];
cx q[12],q[26];
tdg q[26];
cx q[12],q[26];
sdg q[28];
h q[29];
cx q[27],q[29];
tdg q[29];
cx q[13],q[29];
t q[29];
cx q[27],q[29];
cx q[13],q[29];
cx q[13],q[27];
tdg q[27];
cx q[13],q[27];
s q[29];
cx q[26],q[29];
tdg q[29];
cx q[13],q[29];
tdg q[29];
cx q[26],q[29];
cx q[13],q[29];
cx q[13],q[26];
t q[26];
cx q[13],q[26];
h q[29];
cx q[27],q[28];
t q[28];
cx q[12],q[28];
tdg q[28];
cx q[27],q[28];
cx q[12],q[28];
cx q[12],q[27];
t q[27];
cx q[12],q[27];
h q[28];
