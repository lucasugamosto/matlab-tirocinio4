%m-file con il quale vengono assegnati dei valori ai parametri. Questa
%operazione deve essere eseguita prima di simulare il sistema in simulink
%altrimenti i blocchi contenenti le matrici non hanno nessun valore
%definito

%CASO1 : valori con i quali si ottiene il caso di due autovalori immaginari
%puri
m1 = 1;
m2 = 1;
k = 2;
c = 0;

x0 = [0;5;0;0];         %x1: posizione iniziale della massa m1
                        %x2: posizione iniziale della massa m2 (x1+L)
                        %x3: velocità iniziale della massa m1
                        %x4: velocità iniziale della massa m2

A = [0 0 1 0;0 0 0 1;-k/m1 k/m1 -c/m1 c/m1;k/m2 -k/m2 c/m2 -c/m2];
B = [0;0;1/m1;0];
C = [1 0 0 0];
D = 0;
value = 1;

[As,Bs,Cs,Ds] = CreazioneCompensatore1(A,B,C,value);

%funzione che calcola il numeratore e denominatore della funzione di
%trasferimento dell'impianto P(s)
[Np,Dp] = ss2tf(A,B,C,D);

%funzione che calcola il numeratore e denominatore della funzione di
%trasferimento del compensatore C(s)
[Nc,Dc] = ss2tf(As,Bs,Cs,Ds);