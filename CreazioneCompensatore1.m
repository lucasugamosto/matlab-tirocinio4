function CreazioneCompensatore1(A,B,C,value)
    %E' possibile progettare un compensatore che stabilizza asintoticamente
    %il sistema a ciclio chiuso mediante retroazione dinamica dall'uscita
    %se e solo se il sistema S è STABILIZZABILE e RILEVABILE.
    syms s;
    dim_A = size(A);
    n = dim_A(1);
    I = eye(n);
    
    %CONTROLLO STABILIZZABILITA' DEL SISTEMA S
    P = [];
    
    for i = 0:n-1
        P = horzcat(P,(A^i)*B);
    end
    
    if rank(P) == n
        fprintf("il sistema è raggiungibile e quindi stabilizzabile\n")
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
        fprintf("il sistema non è raggiungibile ma è stabilizzabile\n");
        %quindi tutti e soli gli autovalori che fanno parte del
        %sottosistema irraggiungibile sono quelli con parte reale minore di
        %0
    end
    
    %CONTROLLO RILEVABILITA' DEL SISTEMA S
    Q = [];
    
    for i = 0:n-1
        Q = vertcat(Q,C*(A^i));
    end
    
    if rank(Q) == n
        fprintf("il sistema è osservabile e quindi rilevabile\n");
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
        fprintf("il sistema non è osservabile ma è rilevabile\n");
        %quindi tutti e soli gli autovalori che fanno parte del
        %sottosistema inosservabile sono quelli con parte reale minore di
        %0
    end
        
    %PROGETTAZIONE DI Cm(s)
    %In questo caso l'inseguimento di traiettoria è del tipo r(t) = M1*t+M0
    %quindi r(s) ha polo 0 con molteplicità 2. 
    %La funzione di trasferimento P(s) del sistema ha già quel polo con la
    %stessa molteplicità, quindi Cm(s) = 1.
    
    %Siccome Cm(s) = 1 allora le matrici Am,Bm,Cm,Dm sono così definite:
    Am = 0;
    Bm = 0;
    Cm = 0;
    Dm = 1;
    
    %PROGETTAZIONE DI Cs(s)
    fprintf("gli autovalori della matrice A sono:\n");
    autovalori_A = eig(A)
    
    F = CalcoloMatriceF(A,B,value,rag);
    V = CalcoloMatriceV(A,C,value+1,oss);
    
    %le matrici As,Bs,Cs,Ds sono di seguito definite:
    As = A-(V*C)+(B*F);
    Bs = V;
    Cs = -F;
    Ds = 0;
    
    
end