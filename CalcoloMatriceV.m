function V = CalcoloMatriceV(A,C,value,oss)
    %funzione che calcola la matrice V per il sistema S calcolando la
    %corrispettiva matrice F per il sistema duale Sd
    syms t
    dim_A = size(A);
    n = dim_A(1);
    I = eye(n);
    
    if oss == 1
        %CASO sistema osservabile => posso cambiare tutti gli autovalori
        %della matrice A
        
        %calcolo della matrice di raggiungibilità per il sistema duale Sd
        Ad  = A';
        Bd = C';
        
        Pd = [];
        for i = 0:n-1
            Pd = horzcat(Pd,(Ad^i)*Bd);
        end
        
        %calcolo del vettore tau'
        vet = [];
        
        for i = 1:n
            if i == n
                vet = horzcat(vet,1);
            else
                vet = horzcat(vet,0);
            end
        end
        
        tau = vet*(inv(Pd));
        
        %calcolo dei nuovi autovalori per la matrice Ad+BdFd = A-VC
        autoval_Ad = eig(Ad)
        newAutoval = [];
        
        for i = 1:n
            val = real(autoval_Ad(i));
            val = round(val);
            val = val-value;
            newAutoval = vertcat(newAutoval,val);
        end
        newAutoval
        
        %calcolo del polinomio caratteristico desiderato Pdes
        Pdes = 1;
        for i = 1:n
            val = t-newAutoval(i);
            Pdes = Pdes*val;
        end
        Pdes = expand(Pdes);
        vet = sym2poly(Pdes);
        
        %calcolo della matrice Fd
        mat = zeros(n);
        j = 0;
        for i = n+1:-1:1
            val = vet(i)*(A^j);
            mat = mat+val;
            j = j+1;
        end
        
        Fd = -(tau*mat)
        
        %calcolo della matrice V
        V = -(Fd')
    
    else
        %CASO sistema non osservabile ma rilevabile, allora studio degli
        %autovalori della sola parte osservabile
        
        %calcolo delle matrici T e T^-1 per la forma di kalman
        
        %calcolo della matrice di osservabilità Q e la matrice duale Pd
        Q = [];
        for i = 0:n-1
            Q = vertcat(Q,C*(A^i));
        end
        
        Ad = A';
        Bd = C';
        Pd = Q';
        
        %calcolo una base per il sottospazio Xr per il sistema duale Sd
        Xr = [];
        dim_Xr = 0;
        nr = rank(Pd);
        j = 1;
        
        while dim_Xr < nr
            if Pd(:,j) == 0
                j = j+1;
            else
                Xr = horzcat(Xr,Pd(:,j));
                dim_Xr = dim_Xr+1;
                j = j+1;
            end
        end
        
        %calcolo del complemento di Xr
        j = 1;
        while dim_Xr < n
            p = horzcat(Xr,I(:,j));
            if rank(p) == n
                dim_Xr = dim_Xr+1;
                j = j+1;
            else
                p = Xr;
                j = j+1;
            end
        end
        T_inv = p;
        T = inv(T_inv);
        
        %calcolo delle nuove matrici nella forma di kalman
        newAd = T*Ad*T_inv;
        newBd = T*Bd;
        
        %calcolo delle matrici della sola parte raggiungibile
        Adr = [];
        Bdr = [];
        
        for i = 1:nr
            Adr = horzcat(Adr,newAd(1:nr,i));
            Bdr = vertcat(Bdr,newBd(i,:));
        end
        
        %calcolo degli autovalori appartenenti al sottosistema
        %raggiungibile per il sistema duale e dei nuovi autovalori
        autoval_Adr = eig(Adr);
        newAutoval = [];
        
        for i = 1:nr
            val = real(autoval_Adr(i));
            val = round(val);
            val = val-value;
            newAutoval = vertcat(newAutoval,val);
        end
        
        %calcolo della matrice di raggiungibilità per il sottositema
        %raggiungibile del sistema duale Sd
        Pdr = [];
        
        for i = 0:nr-1
            Pdr = horzcat(Pdr,(Adr^i)*Br);
        end
        
        %calcolo del vettore tau' e del polinomio caratteristico desiderato
        %Pdes
        vet = [];
        
        for i = 1:nr
            if i == nr
                vet = horzcat(vet,1);
            else
                vet = horzcat(vet,0);
            end
        end
        
        tau = vet*(inv(Pdr));
        
        Pdes = 1;
        for i = 1:nr
            val = t-newAutoval(i);
            Pdes = Pdes*val;
        end
        Pdes = expand(Pdes);
        vet = sym2poly(Pdes);
        
        mat = zeros(nr);
        j = 0;
        for i = nr+1:-1:1
            val = vet(i)*(Adr^j);
            mat = mat+val;
            j = j+1;
        end
        
        %calcolo della matrice F
        Fdr = -tau*mat;
        vet = zeros(n-nr,1);
        
        while nr < n
            Fdr = horzcat(Fdr,vet);
            nr = nr+1;
        end
        
        newFd = Fdr;
        Fd = newFd*T;
        
        V = -(Fd');
    end 
end