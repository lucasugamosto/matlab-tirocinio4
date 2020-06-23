function [Ncom,Dcom] = CreazioneCompensatore2(A,B,C,D,value,OM)
    %E' possibile progettare un compensatore che stabilizza asintoticamente
    %il sistema a ciclio chiuso mediante retroazione dinamica dall'uscita
    %se e solo se il sistema S è STABILIZZABILE e RILEVABILE.
    
    %progettazione per il caso numerico di due autovalori immaginari puri
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
    %In questo caso l'inseguimento di traiettorie è del tipo
    %r(t) = M*sin(OM*t + f);
    
    %Gli autovalori di A sono: 0,0,2j,-2j => scelgo OM != 2;
    
    %Allora definisco Cm(s) = 1/fi(s) con fi(s) = (s^2 + OM^2).
    if OM == 2
        fprintf("il valore di OMEGA scelto è errato, riprovare\n");
        return
    else
        Am = [-1 0 -(OM^2);1 0 0;0 1 0];
        Bm = [1;0;0];
        Cm = [0 0 1];
        Dm = 0;
    end
    
    %ora si trova una rappresentazione dello spazio di stato per la serie
    %di Cm(s) e P(s)
    dim_Am = size(Am);
    mat = zeros(dim_Am(1),dim_A(2));
    
    Aseg = [Am mat;B*Cm A];
    Bseg = [Bm;B*Dm];
    Cseg = [D*Cm C];
    Dseg = [D*Dm];
    
    %PROGETTAZIONE DI Cs(s)
    %ora si progetta lo stabilizzatore per il sistema (Aseg,Bseg,Cseg,Dseg)
    fprintf("gli autovalori della matrice Aseg sono:\n");
    autovalori_Aseg = eig(Aseg);
    
    F = CalcoloMatriceF(Aseg,Bseg,value,rag);
    V = CalcoloMatriceV(Aseg,Cseg,value+1,oss);
    
    As = Aseg - (V*Cseg) + (Bseg*F);
    Bs = V;
    Cs = -F;
    Ds = 0;
    
    %ora si scrive la rappresentazione nello spazio di stato per C(s) che è
    %la serie di Cs(s) e Cm(s)
    dim_As = size(As);
    mat = zeros(dim_As(1),dim_Am(2));
    
    Ac = [As mat;Bm*Cs Am];
    Bc = [Bs;Bm*Ds];
    Cc = [Dm*Cs Cm];
    Dc = [Dm*Ds];

    [Ncom,Dcom] = ss2tf(Ac,Bc,Cc,Dc);
end