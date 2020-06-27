%m-file con il quale vengono assegnati i valori ai parametri. Questa
%operazione deve essere eseguita prima di simulare il sistema in simulink
%altrimenti i blocchi contenenti le matrici o altro non hanno nessun valore
%definito
fprintf("CASO: 2 autovalori immaginari puri\n");

%caso1 : quaterna di valori m1,m2,k,c tali che gli autovalori della matrice
%A diversi da 0 siano immaginari puri (+-wj)
m1 = 1;
m2 = 1;
k = 2;
c = 0;

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

fprintf("-----------------------------------------------------------\n");
fprintf("PROGETTAZIONE DEL COMPENSATORE NEL CASO DI UNA TRAIETTORIA\n DEL TIPO r(t) = M1t+M0\n");

fprintf("\ngli autovalori della matrice A sono:\n");
autovalori_A = eig(A)

fprintf("\nnumeratore e denominatore della f.d.t del sistema:\n");
[Nsys,Dsys] = ss2tf(A,B,C,D)

[Ncom1,Dcom1] = CreazioneCompensatore1(A,B,C,1);
fprintf("\nnumeratore e denominatore della f.d.t del compensatore:\n");
Ncom1
Dcom1

fprintf("-----------------------------------------------------------\n");
fprintf("PROGETTAZIONE DEL COMPENSATORE NEL CASO DI UNA TRAIETTORIA\n DEL TIPO r(t) = Msin(OMt+fi)\n");

[Aseg,Ncom2,Dcom2] = CreazioneCompensatore2(A,B,C,D,1,3);

fprintf("\ngli autovalori della matrice Aseg sono:\n");
autovalori_Aseg = eig(Aseg)

fprintf("\nnumeratore e denominatore della f.d.t del compensatore:\n");
Ncom2
Dcom2
