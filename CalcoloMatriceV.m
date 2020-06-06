function V = CalcoloMatriceV(A,C,value)
    %funzione che calcola la matrice V risolvendo il problema di trovare la
    %matrice F per il sistema duale Sd
    syms t
    dim_A = size(A);
    n = dim_A(1);
    
    %calcolo delle matrici duali Ad e Bd
    Ad = A';
    Bd = C';
    
    %calcolo della matrice di raggiungibilità Pd
    Pd = [];
    
    for i = 0:n-1
        Pd = horzcat(Pd,(Ad^i)*Bd);
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
    
    tau = vet*(inv(Pd));
    
    %calcolo dei nuovi autovalori
    autovalori_Ad = eig(Ad);
    
    newAutovalori = [];
    for i = 1:n
        val = autovalori_Ad(i)-value;
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
    
    %calcolo della matrice di retroazione Fd
    mat = zeros(n);
    j = 0;
    for i = n+1:-1:1
        val = vet(i)*(Ad^j);
        mat = mat+val;
        j = j+1;
    end

    Fd = -tau*mat;
    
    %calcolo della matrice V
    V = -(Fd');
end