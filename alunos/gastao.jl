#    ÁLGEBRA LINEAR COMPUTACIONAL
#
#    BIBLIOTECA DE FUNÇÕES DO JULIA
#
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------




#    FUNÇÕES BÁSICAS
#-------------------------------------------------------------------------------
function somatorio(v)
    u=ones(length(v))
    return v'u
end
#-------------------------------------------------------------------------------
function media(v)
    return (1/length(v))*somatório(v)
end
#-------------------------------------------------------------------------------
function media_ponderada(v,w)
    if length(v)!=length(w)
        print("A quantidade de dados deve ser a mesma quantidade de pesos.")
    else
        return (1/somatório(w))*v'w
    end
end
#-------------------------------------------------------------------------------
function produto_interno(u,v)
    if length(v)!=length(w)
        print("As quantidades de coordenadas dos dois vetores devem ser as mesmas!!!")
    else
        return u'v
    end
end
#-------------------------------------------------------------------------------
angulo(x,y)=acos(x'*y/(norm(x)*norm(y)))
#-------------------------------------------------------------------------------
function norma(v)
    return sqrt(v'v)
end
#-------------------------------------------------------------------------------
function normalizar(v)
    return v/norma(v)
end
#-------------------------------------------------------------------------------
function distancia(u,v)
    return norma(u-v)
end
#-------------------------------------------------------------------------------
# Encontra, entre os vetores de Y, o vetor que está mais perto de x.
function vetor_mais_perto(x,Y)
    y=Y[argmin([norm(x-y) for y in Y])]
    println("O vetor x é ", x)
    println("O vetor do conjunto Y=", [Y])
    println("mais próximo de x é o vetor ", y)
    return Y[argmin([norm(x-y) for y in Y])]
end

#-------------------------------------------------------------------------------
# Define o i-ésimo vetor canônico de n coordenadas.
function vetor_canonico(i,n)
    if 0<i<=n
        return [zeros(i-1);1;zeros(n-i)]
        else println("Valores inválidos.") 
    end
end
#-------------------------------------------------------------------------------
function combinacao_linear(escalares,vetores)
    return escalares.*vetores
end
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------





#    PROJEÇÃO
#-------------------------------------------------------------------------------
function projeção(a,v)
    v=v/norm(v)
    return (a'v)*v
end
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------




#    DECOMPOSIÇÃO QR
#-------------------------------------------------------------------------------
function decomposição_qr(A)
    n,m=size(A)              
    dim=min(n,m)                # O valor de parada do código é o menor valor entre a quantidade de linhas e colunas de A.
    Q=zeros(n,dim)		
    R=zeros(dim,m)
    for i=1:dim			# Construção das matrizes Q e R.		
        
        Q[:,i]=A[:,i]/norm(A[:,i])	# Coluna i de Q (normalizado)
        
        R[i,:]=Q[:,i]'*A	# Linha i de R: Coeficientes dos vetores de A descritos na base Q por projeção (produto interno)
        
        S=Q[:,i]*R[i,:]'    	# S é a matriz a ser subtraída de A

        A=A-S        
    end
   return Q,R 
end
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------




#    DECOMPOSIÇÃO LU
#-------------------------------------------------------------------------------
function decomposição_lu(A)
    n,m=size(A)
    dim=min(n,m)                # O valor de parada do código é o menor valor entre a quantidade de linhas e colunas de A.
    L=zeros(n,dim)
    U=zeros(dim,m)
    
    for i=1:dim			# Construção das matrizes Q e R.

        L[:,i]=A[:,i]/A[i,i]    # Coluna i de L (Note que o elemento da diagonal principal é 1)

        U[i,:]=A[i,:]	        # Linha i de U
        
        S=L[:,i]*A[i,:]' 	# S é a matriz que iremos subtrair: L[:,i]*U[i,:]

        A=A-S
    end
    
    return L,U
end
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------




#-------------------------------------------------------------------------------
#    MATRIZES TRIANGULARES
#-------------------------------------------------------------------------------
#    SUBSTITUIÇÃO NORMAL
#-------------------------------------------------------------------------------
function retrosubstituicao_triangular_inferior_quadrada(A,b)
    n=length(b)
    x=zeros(n)
    D=0
     
    for i=1:n

        D=A[i,:]'*x
        
        x[i]=(b[i]-D)/A[i,i]
    end
    
    return x 
end
#-------------------------------------------------------------------------------
#    SUBSTITUIÇÃO REVERSA
#-------------------------------------------------------------------------------
function retrosubstituicao_triangular_superior_quadrada(A,b)
    n=length(b)
    x=zeros(n)
    D=0
     
    for i=n:-1:1

        D=A[i,:]'*x
        
        x[i]=(b[i]-D)/A[i,i]
    end
    
    return x 
end
#-------------------------------------------------------------------------------
#    RESOLVER SISTEMA LINEAR POR DECOMPOSIÇÃO LU
#-------------------------------------------------------------------------------
function resolver_LU(A,b)
    
    L,U=decomposição_lu(A)
    y=substituição_triangular_inferior(L,b)
    x=substituição_triangular_superior(U,y)
    
    return x
end
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------





#--------------------------------------------------------------------------------------------------------------------------
#    ESCALONAMENTO COM PIVOTEAMENTO PARCIAL (MATRIZ QUADRADA, PIVO NÃO NULO)
#-------------------------------------------------------------------------------
# Entrada: uma matriz quadrada não singular e um vetor linha b
# Saída: uma matriz quadrada não singular escalonada e um vetor linha b que geram
# um sistema equivalente ao dos dados de entrada.
#-------------------------------------------------------------------------------
function escalonamento_com_pivoteamento(A,b)
    A_amp=[A b']     #Matriz ampliada de A
    n,m=size(A)
  
    for pivo=1:(n-1)
        A_amp=pivoteamento_parcial(A_amp,pivo)
        A_amp=zerar_abaixo(A_amp,pivo)
    end
    
    A=A_amp[:,1:m]
    b=A_amp[:,m+1]
    return A,b
end
#-------------------------------------------------------------------------------
function pivoteamento_parcial(A,pivo)
    Ap=copy(A)
    n_pivo=novo_pivo(Ap[:,pivo],pivo)
    
    linha=Ap[pivo,:]                 #troca as linhas
    Ap[pivo,:]=Ap[n_pivo,:]
    Ap[n_pivo,:]=linha
    
    return Ap
end
#-------------------------------------------------------------------------------
function novo_pivo(coluna,pivo)
    
    maximo=abs(coluna[pivo])
    n_pivo=pivo
        
    for i=pivo:length(coluna)
        if (abs(coluna[i]))>maximo
            maximo=abs(coluna[i])
            n_pivo=i
        end
    end
    
    return n_pivo
end
#-------------------------------------------------------------------------------
function zerar_abaixo(A,pivo)
    n=length(A[:,pivo])

    for i=pivo+1:n
    
        A[i,:]+=-(A[i,pivo]/A[pivo,pivo])A[pivo,:]
    end
    
    return A
end
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------




#-------------------------------------------------------------------------------




#--------------------------------------------------------------------------------------------------------------------------
#    ELIMINAÇÃO GAUSSIANA COM PIVOTEAMENTO PARCIAL (MATRIZ QUADRADA, PIVO NÃO NULO)
#-------------------------------------------------------------------------------
# Preciso das funções: norma; multiplicação matricial; transposta.
# Entrada: uma matriz com os pontos dados
# Saída: o vetor direção "v" da reta que melhor se aproxima dos pontos dados
#-------------------------------------------------------------------------------
function resolver_escalonamento_com_pivoteamento(A,b)
    A,b=escalonamento_com_pivoteamento(A,b)
    x=retrosubstituicao_triangular_superior_quadrada(A,b)
   return x
end
#-------------------------------------------------------------------------------
function retrosubstituicao_triangular_superior_quadrada(A,b)
    n=length(b)
    x=zeros(n)
    D=0
     
    for i=n:-1:1

        D=A[i,:]'*x
        
        x[i]=(b[i]-D)/A[i,i]
    end
    
    return x 
end
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------




#--------------------------------------------------------------------------------------------------------------------------
#    SVD
#-------------------------------------------------------------------------------
# Preciso das funções: norma; multiplicação matricial; transposta.
# Entrada: uma matriz com os pontos dados
# Saída: o vetor direção "v" da reta que melhor se aproxima dos pontos dados
#-------------------------------------------------------------------------------
function SVD(pontos)
    A=copy(pontos) #matriz criada com os pontos dados (Tem que testar se precisa pegar a transposta.)
    v=A[:,1]       
    v=v/norm(v)    #Primeiro candidato a "v" é o primeiro ponto (normalizado) dado (vendo os pontos como os vetores colunas de A).
    erro=0         
    erro_novo=norm(A-v*(A'*v)')  #O erro é verificado vendo a "distância" entre os pontos dados e a reta dada pelo vetor v".
    
    while(abs(erro-erro_novo)>0.01)    # A cada passo o erro diminui mas não necessariarmente converge para zero. 
                                       # Queremos que o processo pare quando a diferença do novo erro para o antigo seja "pequeno".
        erro=erro_novo
        
        for i=1:2
            v=A'*v
            v=v/norm(v)
            A=A'
            end                         # Depois desses dois processos temos o novo candidato a "melhor" v!!!
         
        erro_novo=norm(A-v*(v'*A))      
        
        println(erro)
        println(erro_novo)
        println(erro-erro_novo)         #Deixei os prints pra ver a variação do erro!!!
        println()
    end
    

    return v                            # "Melhor" v!!!    
end
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------




#--------------------------------------------------------------------------------------------------------------------------
#    SISTEMAS DINÂMICOS LINEARES
#------------------------------------------------------------------------------------------
# A função dá o resultado de k iterações da matriz A aplicadas a partir do vetor x0
# Entradas: 1 matriz, 1 vetor, 1 inteiro positivo
# Saída: 1 vetor
#------------------------------------------------------------------------------------------
function dinamica(A,x0,k)
    x=x_0		# Dado inicial.

    for i=1:k		# k iterações.
        x=A*x		# Matriz a sendo aplicada no resultado da iteração anterior.
    end

    return x		# Vetor resultado de todas as iterações.
end
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------




#--------------------------------------------------------------------------------------------------------------------------
#    INTERPOLAÇÃO
#------------------------------------------------------------------------------------------
# Encontrando um candidatoa solução (Necessário para a função "interpolação"!).
function candidato_solucao(pontos,grau,quantidade)
    X=ones(quantidade,grau+1)
        
        for i=1:grau
            X[:,i+1]=X[:,i].*pontos[:,1]
        end
    
    Y=pontos[:,2]

    A=X[1:(grau+1),:]\Y[1:(grau+1)]
    
   return X,A,Y 
end
#------------------------------------------------------------------------------------------
function interpolacao(pontos,grau)
    quantidade=length(pontos[:,1])
    if quantidade<grau+1
        variaveis_livres=grau+1-quantidade
        println("O problema possui infinitas soluções como ", variaveis_livres, " variáveis livres.")
        return variaveis_livres
       
    else
        
        X,A,Y=candidato_solucao(pontos,grau,quantidade)
        
        if quantidade==grau+1
            return A
        
        else
            if norm(X*A-Y)<0.00001
                return A
            else
                return "O problema não tem solução!!!"
            end
        end
    end
end
#------------------------------------------------------------------------------------------
# Polinômio interpolador
a=interpolacao(pontos,grau)
p(x)=a[1]+sum( a[i] * x^(i-1) for i = 2:length(a))
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

