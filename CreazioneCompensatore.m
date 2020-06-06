function [F,V,Ak] = CreazioneCompensatore(A,B,C,value)
    %funzione che calcola la matrice dinamica del compensatore che permette
    %di ottenere stabilit� asintotica del sistema a ciclo chiuso
    
    F = CalcoloMatriceF(A,B,value);
    V = CalcoloMatriceV(A,C,value+2);
    
    %value � il valore che indica di quanto devono variare i nuovi 
    %autovalori rispetto a quelli della matrice dinamica A.
    
    %Gli autovalori di A+BF sono scelti in modo che i corrispondenti modi
    %siano tutti convergenti a zero, ma in maniera pi� rapida rispetto ai
    %modi naturali
    
    %Gli autovalori di A-VC sono scelti in modo che i corrispondenti modi
    %siano tutti convergenti a zero, pi� rapidamente rispetto a quelli di
    %A+BF, per far si che la convergenza a zero dell'errore di stima sia
    %pi� rapida.
    
    Ak = A-(V*C)+(B*F);
end