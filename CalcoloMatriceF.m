function F = CalcoloMatriceF(A,B,value)
    %funzione che permette di calcolare la matrice F utilizzando la formula
    %di Ackermann, favorito rispetto alla formula di Mitter poichè B ha una
    %colonna
    syms t
    dim_A = size(A);
    n = dim_A(1);
    
    %calcolo matrice di raggiungibilità P
    P = [];
    for i = 0:n-1
        P = horzcat(P,(A^i)*B);
    end
        
    %calcolo del vettore tau
    vet = [];
    for i = 1:n
        if i == n
            vet = horzcat(vet,1);
        else
            vet = horzcat(vet,0);
        end
    end
        
    tau = vet*inv(P);
        
    %calcolo dei nuovi autovalori
    autovalori_A = eig(A);
        
    newAutovalori = [];
    for i = 1:n
        val = autovalori_A(i)-value;
        newAutovalori = vertcat(newAutovalori,val);
    end
        
    %calcolo del polinomio caratteristico desiderato Pdes
    Pdes = 1;
    for i = 1:n
        val = t-newAutovalori(i);
        Pdes = Pdes*val;
    end
    Pdes = expand(Pdes);
    vet = sym2poly(Pdes);
    
    %calcolo della matrice di retroazione F
    mat = zeros(n);
    j = 0;
    for i = n+1:-1:1
        val = vet(i)*(A^j);
        mat = mat+val;
        j = j+1;
    end

    F = -tau*mat;
end