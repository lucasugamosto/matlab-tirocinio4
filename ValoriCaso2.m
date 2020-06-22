%m-file con il quale vengono assegnati i valori ai parametri. Questa
%operazione deve essere eseguita prima di simulare il sistema in simulink
%altrimenti i blocchi contenenti le matrici o altro non hanno nessun valore
%definito

%caso2 : quaterna di valori m1,m2,k,c tali che gli autovalori della matrice
%A diversi da 0 siano complessi coniugati
m1 = 2;
m2 = 2;
k = 10;
c = 6;

x = [0;10;0;0];
%x è il vettore contenente lo stato iniziale del sistema, in esso:
%x1 = posizione iniziale della massa m1;
%x2 = posizione iniziale della massa m2 che si trova a distanza L da m1;
%x3 = velocità iniziale della massa m1;
%x4 = velocità iniziale della massa m1;

A = [0 0 1 0;0 0 0 1;-k/m1 k/m1 -c/m1 c/m1;k/m2 -k/m2 c/m2 -c/m2];
B = [0;0;1/m1;0];
C = [1 0 0 0];
D = 0;

autovalori_A = eig(A);

syms s
dim_A = size(A);
n = dim_A(1);
I = eye(n);

P = D + C*(inv(s*I-A))*B;
pretty(P)

[Ncom,Dcom] = CreazioneCompensatore1(A,B,C,1);