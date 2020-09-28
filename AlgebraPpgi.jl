#------------------------------------------------------------------------------------------
# Retorna o maior valor de um vetor
# Entradas: 1 vetor
# Saídas: 1 valor numérico
#------------------------------------------------------------------------------------------
function valor_maximo(v)
    m = -Inf
    for i in v
        if(i > m) m = i end
    end
    return m
end


#------------------------------------------------------------------------------------------
# Garante que o retorno será um vetor caso a entrada seja uma matriz de 1 coluna/linha
# Entradas: 1 vetor/matriz
# Saídas: 1 vetor
#------------------------------------------------------------------------------------------
function garante_vetor(v) 
    @assert ndims(v) == 1 || (valor_maximo(size(v)) == length(v))
    return v[:]
end


#------------------------------------------------------------------------------------------
# Garante que o retorno será uma matriz de 1 coluna caso a entrada seja um vetor
# Se v tiver mais de 2 dimensões, causará um erro
# Entradas: 1 vetor
# Saídas: 1 matriz
#------------------------------------------------------------------------------------------
function garante_matriz(v) 
    @assert 1 <= ndims(v) <= 2 
    return v[:,:]
end


#------------------------------------------------------------------------------------------
# Soma de duas matrizes
# Entradas: 2 matrizes
# Saídas: 1 matriz
#------------------------------------------------------------------------------------------
function soma_matrizes(A, B)
    A = garante_matriz(A)
    B = garante_matriz(B)

    @assert size(A) == size(B) && ndims(A) == 2

    mA, nA = size(A)
    result = zeros(mA, nA)

    for i = 1:mA
        for j = 1:nA
            result[i, j] = A[i, j] + B[i, j]
        end
    end
    return result
end


#------------------------------------------------------------------------------------------
# Multiplicação de duas matrizes
# Entradas: 2 matrizes
# Saídas: 1 matriz
#------------------------------------------------------------------------------------------
function multiplicar_matrizes(A, B)
    A = garante_matriz(A)
    B = garante_matriz(B)

    mA, nA = size(A)
    mB, nB = size(B)
        
    @assert mA == nB

    result = zeros(mA, nB)

    for i = 1:mA
        for j = 1:nB
            for k = 1:nA
                result[i, j] += A[i, k] * B[k, j]
            end
        end
    end
    return result
end


#------------------------------------------------------------------------------------------
# Produto de matrizes (ou recursiva) recebe um número qualquer de matrizes e 
# retorna o produto entre estas, calculado de forma recursiva
# Entradas: "n" matrizes
# Saídas: 1 matriz
#------------------------------------------------------------------------------------------
function prod_mat_rec(A, B...) # B... significa um número variável de argumentos
    if(length(B) == 0) # se não houver nenhum argumento adicional, retorna a própria matriz A
        return A
    else # se houver algum argumento adicional, retorna a matriz A multiplicada pelo produto das matrizes restantes
        return multiplicar_matrizes(A, prod_mat_rec(B...))
    end
end


#------------------------------------------------------------------------------------------
# Soma de matrizes (ou recursiva) recebe um número qualquer de matrizes e 
# retorna a soma entre estas, calculado de forma recursiva
# Entradas: "n" matrizes
# Saídas: 1 matriz
#------------------------------------------------------------------------------------------
function soma_mat_rec(A, B...) # B... significa um número variável de argumentos
    if(length(B) == 0) # se não houver nenhum argumento adicional, retorna a própria matriz A
        return A
    else # se houver algum argumento adicional, retorna a matriz A multiplicada pelo produto das matrizes restantes
        return mat_sum(A, soma_mat_rec(B...))
    end
end


#------------------------------------------------------------------------------------------
# Produto interno entre dois vetores
# Entradas: 2 vetores
# Saídas: 1 valor numérico
#------------------------------------------------------------------------------------------
function produto_interno(u, v)
    u1 = garante_vetor(u)
    v1 = garante_vetor(v)

    result = 0

    for i = 1:length(u)
        result += u[i] * v[i]
    end

    return result
end


#------------------------------------------------------------------------------------------
# Retorna a norma de um vetor (comprimento)
# Entradas: 1 vetores
# Saídas: 1 valor numérico
#------------------------------------------------------------------------------------------
function norma(a) 
    return sqrt(produto_interno(a, a))
end 


#------------------------------------------------------------------------------------------
# Retorna a norma do máximo valor de um vetor
# Entradas: 1 vetore
# Saídas: 1 valor numérico
#------------------------------------------------------------------------------------------
function norma_do_maximo(a) 
    return valor_maximo(abs.(a))
end 


#------------------------------------------------------------------------------------------
# Testa a ortogonalidade entre dois vetores
# Entradas: 2 vetores
# Saídas: booleano (true ou false)
#------------------------------------------------------------------------------------------
function dois_vetores_ortogonais(a, b, tolerancia=10e-5)
    return produto_interno(a, b) <= tolerancia
end


#------------------------------------------------------------------------------------------
# Testa a ortogonalidade entre um vetor e o plano
# Entradas: 1 vetor, 1 matriz que constrói o plano
# Saídas: booleano (true ou false)
#------------------------------------------------------------------------------------------
function vetor_ortogonal_plano(v, P, tolerancia=10e-5)
    return norma(multiplicar_matrizes(P, v)) <= tolerancia
end


#------------------------------------------------------------------------------------------
# Testa se duas retas são ortogonais
# Entradas: 2 retas. Cada reta é representada por 2 vetores.
# "a" e "b" são tuplas com 2 vetores cada uma, que representam o início e o final de segmentos de retas
# Saídas: booleano (true ou false)
#------------------------------------------------------------------------------------------
function retas_ortogonais(a, b)
    x1, x2 = a; # x1 receberá o vetor do início do segmento da reta "a", e x2 receberá o vetor do final
    y1, y2 = b; # y1 receberá o vetor do início do segmento da reta "b", e y2 receberá o vetor do final
    return dois_vetores_ortogonais(x2 - x1, y2 - y1)
end


#------------------------------------------------------------------------------------------
# Calcula a projeção de um vetor em uma reta
# Entradas: um vetor "v" e uma reta "R" (matriz)
# Saídas: float
#------------------------------------------------------------------------------------------
function projecao_vetor_reta(v, R)
    x1, x2 = R; # x1 receberá o vetor do início do segmento da reta "R", e x2 receberá o vetor do final
    return produto_interno(v, x2 - x1) / (norma(x2 - x1) ^ 2) * (x2 - x1);
end


#------------------------------------------------------------------------------------------
# Calcula a projeção de um vetor em um plano
# Entradas: 2 vetores, 1 é o vetor a ser projetado e o segundo é ortogonal ao plano
# Saídas: float
#------------------------------------------------------------------------------------------
function projecao_vetor_plano(v, u)
    return v - produto_interno(v, u) / (norma(u) ^ 2) * u;
end


#------------------------------------------------------------------------------------------
# Realiza a eliminação gaussiana, usa uma sequência de operações de linha elementares para 
# modificar a matriz até que o canto esquerdo inferior da matriz seja preenchido com zeros, 
# tanto quanto possível
# 
# Entradas: 1 matriz
# Saídas: 1 matriz
#------------------------------------------------------------------------------------------
function eliminacao_gaussiana(A)
    rows = size(A,1)
    cols = size(A,2)
    fraction = A[1,1]/A[1,1]
    # Row index
    row = 1
    U = copy!(similar(A, typeof(fraction)), A)
    # Main loop going through all columns
    for col = 1:(cols-1)
        
        # finding the maximum element for each column
        max_index = argmax(abs.(A[row:end,col])) + row-1
        # Check to make sure matrix is good!
        if (U[max_index, col] == 0)
            continue
        end

        # swap row with highest value for that column to the top
        temp_vector = U[max_index, :]
        U[max_index, :] = U[row, :]
        U[row, :] = temp_vector
        
        # Loop for all remaining rows
        for i = (row+1):rows
            # finding fraction
            fraction = U[i,col]/U[row,col]
            # loop through all columns for that row
            for j = (col+1):cols
                 # re-evaluate each element  
                 U[i,j] -= U[row,j]*fraction
            end
            # Set lower elements to 0
            U[i,col] = 0
        end
        row += 1
    end
    return U
end


#------------------------------------------------------------------------------------------
# Verifica se dois vetores são perpendiculares
# 
# Entradas: 2 vetores
# Saídas: booleano
#------------------------------------------------------------------------------------------
function vetores_perpendiculares(u,v)
    result = produto_interno(u,v)
    if (result == 0) 
        return true
    else
        return false
    end
end


#------------------------------------------------------------------------------------------
# Calcula o ângulo entre dois vetores
# 
# Entradas: 2 vetores
# Saídas: float
#------------------------------------------------------------------------------------------
function angulo_vetores(u,v)
    u1 = garante_vetor(u)
    v1 = garante_vetor(v)
    result = 0
    result = produto_interno(u,v)/(norma(u1)*norma(v1))
    return result
end


#------------------------------------------------------------------------------------------
# Decomposição QR de uma matriz A.
#
# Entrada: uma matriz A.
# Saída: uma matriz ortogonal Q e uma matriz R.
# Autor: Gastão.
#------------------------------------------------------------------------------------------
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


#------------------------------------------------------------------------------------------
# Escalonamento com pivoteamento parcial de uma matriz quadrada não singular com pivôs não-nulos.
#
# Entrada: uma matriz quadrada não singular e um vetor linha b.
# Saída: uma matriz quadrada não singular escalonada e um vetor linha b que geram
# um sistema equivalente ao dos dados de entrada.
# Autor: Gastão.
#------------------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------------------
function pivoteamento_parcial(A,pivo)
    Ap=copy(A)
    n_pivo=novo_pivo(Ap[:,pivo],pivo)
    
    linha=Ap[pivo,:]                 #troca as linhas
    Ap[pivo,:]=Ap[n_pivo,:]
    Ap[n_pivo,:]=linha
    
    return Ap
end
#------------------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------------------
function zerar_abaixo(A,pivo)
    n=length(A[:,pivo])

    for i=pivo+1:n
    
        A[i,:]+=-(A[i,pivo]/A[pivo,pivo])A[pivo,:]
    end
    
    return A
end


#------------------------------------------------------------------------------------------
# Substituição reversa em uma matriz triangular superior quadrada não-singular.
#
# Entrada: uma matriz quadrada triangular não singular A e um vetor linha b
# Saída: o vetor x, solução de Ax=b.
# Autor: Gastão.
#------------------------------------------------------------------------------------------
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


#------------------------------------------------------------------------------------------
# Eliminação gaussiana com pivoteamento parcial de uma matriz quadrada não singular com pivôs não-nulos.
#
# Entrada: uma matriz com os pontos dados.
# Saída: o vetor direção "v" da reta que melhor se aproxima dos pontos dados.
# Autor: Gastão.
#------------------------------------------------------------------------------------------
function resolver_escalonamento_com_pivoteamento(A,b)
    A,b=escalonamento_com_pivoteamento(A,b)
    x=retrosubstituicao_triangular_superior_quadrada(A,b)
   return x
end
#------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------
# Calcula o posto de uma Matriz, que corresponte ao número de linhas ou colunas 
# linearmente independentes da matriz.
# 
# Entradas: 1 vetores.
# Saídas: inteiro.
#------------------------------------------------------------------------------------------
function posto_matriz(A)
    U = eliminacao_gaussiana(A)
    rows = size(U,1)
    rank = rows
    for i = 1:rows
        row_sum =0
        for j = 1:length(U[i, :])
            row_sum += U[i,j]
        end
        if (row_sum ==0) 
            rank-=1
        end
    end
    return rank
end


#------------------------------------------------------------------------------------------
# SVD: A função dá o vetor direção da "melhor" reta que se aproxima dos pontos dados.
#
# Entrada: uma matriz com os pontos dados
# Saída: o vetor direção "v" da reta que melhor se aproxima dos pontos dados
# Pré-funções: norma.
# Autor: Gastão.
#------------------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------------------
function teste_SVD()      # Teste básico ainda!!!
    pontos=[1 2 3; 2 4 6 ; 3 6 9; 4 8 12]
    v=SVD(pontos')
    println("O vetor v é ", v)
    a=v=v/v[1]
    println(a)
    println()

    pontos=[1 2 3; 2 4 6 ; 3 6 9; 4 8 15]
    v=SVD(pontos')
    println("O vetor v é ", v)
    a=v=v/v[1]
    println(a)
    println()
end


#------------------------------------------------------------------------------------------
# SISTEMAS DINÂMICOS LINEARES: A função dá o resultado de k iterações da matriz A aplicadas 
# a partir do vetor x0.
#
# Entradas: 1 matriz quadrada, 1 vetor, 1 inteiro positivo.
# Saída: 1 vetor.
# Pré-funções: norma.
# Autor: Gastão.
#------------------------------------------------------------------------------------------
function dinamica(A,x0,k)
    x=x0               # Dado inicial.

    for i=1:k           # k iterações.
        x=A*x           # Matriz a sendo aplicada no resultado da iteração anterior.
    end

    return x            # Vetor resultado de todas as iterações.
end
#------------------------------------------------------------------------------------------
function teste_dinamica(k)
    Tudo_certo=true

    # Uma dinâmica aplicada no vetor nulo sempre dará como resposta o vetor nulo.
    for n=2:k
        A=randn(n,n)                # Matriz A qualquer
        x0=zeros(n,1)               # Dado inicial nulo
        for i=1:k
            x=dinamica(A,x0,i)      # Resultado da dinâmica
            if norm(x)>0.00001      # Verificação se o resultado continua sendo o vetor nulo
                Tudo_certo=false    
            end
        end
    end

    # Se a matriz for a identidade então ela não irá alterar o vetor.
    for n=2:k
        A=zeros(n,n)                  
        for i=1:n A[i,i]=1 end       # Matriz A identidade
        x0=randn(n,1)                # Dado inicial qualquer
        for i=1:k
            x=dinamica(A,x0,i)       # Resultado da dinâmica
            if norm(x0-x)>0.00001    # Verifica se o resultado continua sendo o dado inicial
                Tudo_certo=false
            end
        end
    end

    return Tudo_certo   
end

#-------------------------------------------------------------------------------
#   TRANSPOSTA
#
# Recebe: A -> matriz ou vetor que sofererá transposição
# Retorna transposta -> matriz ou vetor transposto de A
#-------------------------------------------------------------------------------

function transposta(A)
    
    if size(A) == (1,) #verifica se A é um vetor unitário. Se sim, retorna ele mesmo,
        return(A)      #pois não há transposição a ser efetuada
    else
        m,n = size(A)  #calculando o tamanho da matriz/vetor A
        transposta = zeros(n,m) #Gerando uma matriz/vetor com dimensões transpostas preenchido com zeros

        for i=1:m, j=1:n #Para cada coordenada da matriz/vetor A
            transposta[j,i,:] = A[i,j, :] #Substitui o valor n,m da matriz de zeros por m,n da matriz original 
        end
        return transposta# retorna a matriz transposta
    end
end

#-------------------------------------------------------------------------------------
# TESTE TRANSPOSTA
#
# Verifica se matrizes aleatórias de dimensões diferentes são transpostas corretamente.
# A verificação compara cada linha da matriz original, com a coluna da matriz transposta.
# Se a transposição aconteceu corretamente, elas devem ser iguais
#--------------------------------------------------------------------------------------

function verificacao_transposta(A, A_transposta)
    
    if size(A) == (1,)             #verifica se A é um vetor unitário.
        return A == a_transposta   #Se sim, o compara diretamente com a transposta
    else
    
        m,n = size(A)

        for i=1:m  #para cada coluna de A
            if A[i,:] == A_transposta[:,i] #Se a coluna de A é igual à linha de A_transposta

            else
                return false #Senão, não é a transposta
                break
            end
        end
        return true #Caso todas as colunas de A sejam iguais às linhas de A_transposta
    end
end

#-------------------------------------------------------------------------------
#   INVERSA
#
# Recebe: A -> matriz ou vetor que será invertida
# Retorna inversa -> matriz inversa de A (A^⁻1)
#-------------------------------------------------------------------------------
function inversa(A)

    if size(A) != (1,) #Verificando se a matriz é unitária, nesse caso não há inversa 

        m,n = size(A) #Verificando dimensões da matriz A
        
        if m == n  #Verificando se a matriz é quadrada

            matriz_inversa = zeros(m,m)  #Iniciando uma matriz com as mesmas dimensões de A

            for i=1:m #Para cada coluna de A

                vetor_matriz_identidade = zeros(m,1)
                vetor_matriz_identidade[i,1] = 1 #Vetor da matriz identidade correspondente à coluna de A
                vetor_inversa = resolver_LU(A,vetor_matriz_identidade) #Resolvendo o sistema linear para criar o vetor da inversa
                matriz_inversa[:,i]=vetor_inversa #criando a matriz inversa vetor a vetor

            end
            return matriz_inversa

        else  #Se a matriz não é quadrada, ela não é inversível
            return "Matriz não inversível"
        end
    else
        return "Matriz não inversível"
    end
end


#-------------------------------------------------------------------------------------
# TESTE INVERSA
#
# Verifica se A*A^-1 = I
#--------------------------------------------------------------------------------------
function testa_matriz_inversa(A, inversa_A)
    if norm((A*inversa_A) - I) < 0.000001 # Se a norma de A*A^-1 - I é proxima de 0
        return true
    else
        return false
    end
end