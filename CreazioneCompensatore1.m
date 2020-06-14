function [As,Bs,Cs,Ds] = CreazioneCompensatore1(A,B,C,value)
    %E' possibile progettare un compensatore che stabilizza asintoticamente
    %il sistema a ciclio chiuso mediante retroazione dinamica dall'uscita
    %se e solo se il sistema S è STABILIZZABILE e RILEVABILE.
    syms s;
    dim_A = size(A);
    n = dim_A(1);
    I = eye(n);
    
    %controllo stabilizzabilità del sistema S
    P = [];
    
    for i = 0:n-1
        P = horzcat(P,(A^i)*B);
    end
    
    if rank(P) == n
        fprintf("il sistema è raggiungibile\n");
        rag = 1;
    else
        rag = 0;
        autovalori_A = eig(A);
        
        for i = 1:n
            if real(autovalori_A(i)) >= 0
                fprintf("PBH test dell'autovalore %d\n",autovalori_A(i));
                mat = horzcat((A-autovalori_A(i)*I),B);
                if rank(mat) ~= n
                    fprintf("il sistema non è stabilizzabile ->\n");
                    fprintf("-> impossibile progettare un compensatore\n");
                    return
                end
            end
        end
        fprintf("il sistema è stabilizzabile\n");
    end
    
    %controllo rilevabilità del sistema S
    Q = [];
    
    for i = 0:n-1
        Q = vertcat(Q,C*(A^i));
    end
    
    if rank(Q) == n
        fprintf("il sistema è osservabile\n");
        oss = 1;
    else
        oss = 0;
        autovalori_A = eig(A);
        
        for i = 1:n
            if real(autovalori_A(i)) >= 0
                fprintf("PBH test dell'autovalore %d\n",autovalori_A(i));
                mat = vertcat((A-autovalori_A(i)*I),C);
                if rank(mat) ~= n
                    fprintf("il sistema non è rilevabile ->\n");
                    fprintf("-> impossibile progettare un compensatore\n");
                    return
                end
            end
        end
        fprintf("il sistema è rilevabile\n");
    end
    
    %progettazione di Cm(s)
    %In questo caso l'inseguimento di traiettoria è del tipo r(t) = M1*t+M0
    %quindi r(s) ha polo 0 con molteplicità 2. 
    %La funzione di trasferimento P(s) del sistema ha già quel polo con la
    %stessa molteplicità, quindi Cm(s) = 1.
    
    %progettazione dello stabilizzatore Cs(s)
    F = CalcoloMatriceF(A,B,value,rag);
    V = CalcoloMatriceV(A,C,value+1,oss);
    
    %calcolo le matrici As,Bs,Cs,Ds per lo stabilizzatore
    As = A-(V*C)+(B*F);
    Bs = V;
    Cs = -F;
    Ds = 0;
    
    %siccome Cm(s) = 1 allora le matrici Am,Bm,Cm,Dm valgono come segue:
    %Am = 1;
    %Bm = 1;
    %Cm = 1;
    %Dm = 0;
    
    %la rappresentazione nello spazio di stato per C(s) che è la cascata di
    %Cs(s) e Cm(s) è definita dalle seguenti matrici:
    %dim_As = size(As);
    %n = dim_As(1);
    %z = zeros(n,1);
    %mat = Bm*Cs;
    
    %Ac = [As z;mat Am];
    
    %mat = Bm*Ds;
    
    %Bc = [Bs;mat];
    
    %mat = Dm*Cs;
    
    %Cc = [mat Cm];
    
    %Dc = Dm*Ds;
end