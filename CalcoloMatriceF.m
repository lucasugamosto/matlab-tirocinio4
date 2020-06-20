function F = CalcoloMatriceF(A,B,value,rag)
    %funzione che calcola la matrice di retroazione dello stato F
    syms t
    dim_A = size(A);
    n = dim_A(1);
    I = eye(n);
    
    if rag == 1
        %CASO sistema raggiungibile => posso cambiare tutti gli autovalori
        %della matrice A
        
        %calcolo matrice di raggiungibilità P
        P = [];  
        for i = 0:n-1
            P = horzcat(P,(A^i)*B);
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
        tau = vet*(inv(P));
        
        %calcolo dei nuovi autovalori della matrice A+BF
        autoval_A = eig(A);
        newAutoval = [];
        
        for i = 1:n
            if real(autoval_A(i)) >= 0

                val = real(autoval_A(i));
                val = round(val);
                
                while val > -1
                    val = val-value;
                end
                  
                newAutoval = vertcat(newAutoval,val);
            else
                newAutoval = vertcat(newAutoval,autoval_A(i));
            end
        end
        fprintf("i nuovi autovalori per la matrice A+BF sono:\n");
        newAutoval
        
        %calcolo del polinomio caratteristico desiderato Pdes
        Pdes = 1;
        for i = 1:n
            val = t-newAutoval(i);
            Pdes = Pdes*val;
        end
        Pdes = expand(Pdes);
        vet = sym2poly(Pdes);
        
        %calcolo della matrice F
        mat = zeros(n);
        j = 0;
        for i = n+1:-1:1
            val = vet(i)*(A^j);
            mat = mat+val;
            j = j+1;
        end
        
        F = -(tau*mat);
        
    else
        %CASO sistema non raggiungibile ma stabilizzabile poichè gli
        %autovalori della parte non raggiungibile sono con parte reale
        %negativa
        
        %calcolo delle matrici T e T^-1 per la forma di kalman
        
        %calcolo matrice di raggiungibilità P
        P = [];  
        for i = 0:n-1
            P = horzcat(P,(A^i)*B);
        end
        
        %calcolo una base per il sottospazio Xr
        Xr = [];
        dim_Xr = 0;
        nr = rank(P);
        j = 1;
        
        while dim_Xr < nr
            if P(:,j) == 0
                j = j+1;
            else
                Xr = horzcat(Xr,P(:,j));
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
        newA = T*A*T_inv;
        newB = T*B;
        
        %calcolo delle matrici della sola parte raggiungibile
        Arr = [];
        Br = [];
        
        for i = 1:nr
            Arr = horzcat(Arr,newA(1:nr,i));
            Br = vertcat(Br,newB(i,:));
        end
        
        %calcolo degli autovalori appartenenti al sottosistema
        %raggiungibile e dei nuovi autovalori
        autoval_Arr = eig(Arr);
        newAutoval = [];
        
        for i = 1:nr
            if real(autoval_Arr(i)) >= 0

                val = real(autoval_Arr(i));
                val = round(val);
                
                while val > -1
                    val = val-value;
                end
                  
                newAutoval = vertcat(newAutoval,val);
            else
                newAutoval = vertcat(newAutoval,autoval_Arr(i));
            end
        end
        fprintf("i nuovi autovalori per la matrice Arr+BrFr sono:\n");
        newAutoval
        
        %calcolo della matrice di raggiungibilità per il sottosistema
        %raggiungibile
        Pr = [];
        
        for i = 0:nr-1
            Pr = horzcat(Pr,(Arr^i)*Br);
        end
    
        %calcolo della matrice Fr per il sottosistema raggiungibile
        vet = [];
        
        for i = 1:nr
            if i == nr
                vet = horzcat(vet,1);
            else
                vet = horzcat(vet,0);
            end
        end
        
        tau = vet*(inv(Pr));
        
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
            val = vet(i)*(Arr^j);
            mat = mat+val;
            j = j+1;
        end
        
        %calcolo della matrice F
        Fr = -tau*mat;
        vet = zeros(n-nr,1);
        
        while nr < n
            Fr = horzcat(Fr,vet);
            nr = nr+1;
        end
        
        F = Fr*T;
    end      
end