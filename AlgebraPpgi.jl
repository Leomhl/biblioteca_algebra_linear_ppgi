using LinearAlgebra


#******************************************************************************************
#matrizes para testes
notTriangularMatrix = [-1 1/4 1/4 0;
                        1/4 -1 0 1/4;
                        1/4 0 -1 1/4;
                        0 1/4 1/4 -1]

upperTriangularMatrix = [-1 1/4 1/4 0;
                        0 -1 -1 1/4;
                        0 0 1/4 1/4;
                        0 0 0 -1]

lowerTriangularMatrix = [-1 0 0 0;
                        1/4 -1 0 0;
                        1/4 0 -1 0;
                        0 1/4 1/4 -1]

#******************************************************************************************
# VALOR MAXIMO
#Retorna o maior valor de um vetor
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


#******************************************************************************************
# GARANTE VETOR
# Garante que o retorno será um vetor caso a entrada seja uma matriz de 1 coluna/linha
# Entradas: 1 vetor/matriz
# Saídas: 1 vetor
#------------------------------------------------------------------------------------------
function garante_vetor(v) 
    @assert ndims(v) == 1 || (valor_maximo(size(v)) == length(v))
    return v[:]
end


#******************************************************************************************
# GARANTE MATRIZ 
#Garante que o retorno será uma matriz de 1 coluna caso a entrada seja um vetor
# Se v tiver mais de 2 dimensões, causará um erro
# Entradas: 1 vetor
# Saídas: 1 matriz
#------------------------------------------------------------------------------------------
function garante_matriz(v) 
    @assert 1 <= ndims(v) <= 2 
    return v[:,:]
end


#******************************************************************************************
# SOMA MATRIZ
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


#******************************************************************************************
# SOMA MATRIZ
#Soma de matrizes (ou recursiva) recebe um número qualquer de matrizes e 
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




#******************************************************************************************
# MULTIPLICACAO Multiplicação de duas matrizes
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

#******************************************************************************************
# MULTIPLICACAO Produto de matrizes (ou recursiva) recebe um número qualquer de matrizes e 
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




#******************************************************************************************
# PRODUTO INTERNO Produto interno entre dois vetores
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

#******************************************************************************************
# NORMA Retorna a norma de um vetor (comprimento)
# Entradas: 1 vetores
# Saídas: 1 valor numérico
#------------------------------------------------------------------------------------------
function norma(a) 
    return sqrt(produto_interno(a, a))
end 

#******************************************************************************************
#NORMA VETOR
function NormVector(V)
   n = length(V);
   S = 0;
   for i = 1:n
       S = S + (V[i]^2);
   end
   return sqrt(S);
end
#------------------------------------------------------------------------------------------

#******************************************************************************************
#NORMA MATRIZ
function NormMatrix(A)
   m,n = size(A);
   S = 0;
   for i = 1:m
       for j = 1:n
           S = S + (A[i,j]^2);
       end
   end
   return sqrt(S);
end

#******************************************************************************************
# NORMA DO MAXIMO Retorna a norma do máximo valor de um vetor
# Entradas: 1 vetore
# Saídas: 1 valor numérico
#------------------------------------------------------------------------------------------
function norma_do_maximo(a) 
    return valor_maximo(abs.(a))
end 



#******************************************************************************************
#  ANGULO Calcula o ângulo entre dois vetores
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


#******************************************************************************************
# PERPENDICULAR Verifica se dois vetores são perpendiculares
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


#******************************************************************************************
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




#******************************************************************************************
#TRANSPOSTA
#matriz tal que o elemento a_ij é igual ao elemento a_ji da matriz original
#necessária para verificar ortogonalidade

#Entrada: matriz (qualquer tamanho)
#Etapa 1: verificar se é matriz (senão retorna erro)
#Etapa 2: construir matriz transposta vazia com dimensoes corretas
#Etapa 3: popula matriz transposta onde o novo elemento a_ij é igual ao elemento a_ji da matriz original
#Saída: matriz (transposta)
#Autora: gabriella radke
#-------------------------------------------------------------------------------------

function transposta(matrix)
    #checa se input é matriz
    if isa(matrix, Array{Float64,2}) == true  || isa(matrix, Array{Int64,2}) == true
        #checa se input é matriz quadrada
        numberOfLines, numberOfColumns = size(matrix)
        if numberOfLines == numberOfColumns
            #cria matriz transposta vazia com dimensoes corretas para popularmos depois
            transposta = zeros(numberOfColumns, numberOfLines)
            #o elemento a_ij da nova matriz sera o a_ji da original
            for line in 1:numberOfColumns 
                for column in 1:numberOfLines
                    transposta[line,column] = matrix[column,line]
                end
            end    
        else
            #retorna erro se matriz nao e quadrada
            throw(ArgumentError(" 'matrix' parameter should be a SQUARE matrix"))
        end
    else
        #retorna erro se input não é array
        throw(ArgumentError("'matrix' parameter should be a matrix"))
    end
    return transposta 
end

#teste
#eh_triangular(transposta(upperTriangularMatrix))
transposta(upperTriangularMatrix)




#******************************************************************************************
# RETA ORTOGONAL Testa se duas retas são ortogonais
# Entradas: 2 retas. Cada reta é representada por 2 vetores.
# "a" e "b" são tuplas com 2 vetores cada uma, que representam o início e o final de segmentos de retas
# Saídas: booleano (true ou false)
#------------------------------------------------------------------------------------------
function retas_ortogonais(a, b)
    x1, x2 = a; # x1 receberá o vetor do início do segmento da reta "a", e x2 receberá o vetor do final
    y1, y2 = b; # y1 receberá o vetor do início do segmento da reta "b", e y2 receberá o vetor do final
    return dois_vetores_ortogonais(x2 - x1, y2 - y1)
end


#******************************************************************************************
# VETOR ORTOGONAL Testa a ortogonalidade entre dois vetores
# Entradas: 2 vetores
# Saídas: booleano (true ou false)
#------------------------------------------------------------------------------------------
function dois_vetores_ortogonais(a, b, tolerancia=10e-5)
    return produto_interno(a, b) <= tolerancia
end


#******************************************************************************************
# PLANO VETOR ORTOGONAL Testa a ortogonalidade entre um vetor e o plano
# Entradas: 1 vetor, 1 matriz que constrói o plano
# Saídas: booleano (true ou false)
#------------------------------------------------------------------------------------------
function vetor_ortogonal_plano(v, P, tolerancia=10e-5)
    return norma(multiplicar_matrizes(P, v)) <= tolerancia
end


#******************************************************************************************
#ORTOGONAL MATRIZ - verificar se uma matriz é ortogonal
#Note que ser ortogonal é uma propriedade de matrizes quadradas
#portanto, a matriz do input deverá ser quadrada
#Uma matriz é ortogonal quando suas linhas e colunas são vetores ortogonais
#outra forma de saber é se a transposta vezes a matriz é igual identidade
#iremos utilizar a funcao prod_mat_rec para o produto das matrizes

#Entrada: matriz
#Etapa 1: verifica se entrada é matriz quadrada, senao retorna erro
#Etapa 2: calcula transposta da matriz
#Etapa 3: multiplica transposta pela matriz original e verifica se deu identidade
#Saída: true(é ortogonal) ou false(não é ortogonal)
#Autora: gabriella radke
#-------------------------------------------------------------------------------------

function ortogonal(matrix)
     #checa se input é matriz
    if isa(matrix, Array{Float64,2}) == true  || isa(matrix, Array{Int64,2}) == true
        #checa se input é matriz quadrada
        numberOfLines, numberOfColumns = size(matrix)
        if numberOfLines == numberOfColumns
            matriz_transposta = transposta(matrix)
            #se matrix^T * matriz = Id, é ortogonal
            if matriz_transposta*matrix == I
                return true
            else
                return false
            end
        else
            #retorna erro se matriz nao e quadrada
            throw(ArgumentError(" 'matrix' parameter should be a SQUARE matrix"))
        end
    else
        #retorna erro se input não é array
        throw(ArgumentError("'matrix' parameter should be a matrix"))
    end
end


#teste
ortogonal([1 0; 0 1])


#******************************************************************************************
#ORTONORMAL BASE - verificar se um dado conjunto de vetores forma base ortonormal
#um conjunto de vetores forma base ortonormal se cada par de vetores é ortogonal entre si
#e o comprimento de cada vetor é 1
#como existe a funcao dois_vetores_ortogonais que avalia se dois vetores sao ortogonais
#ja implementada, iremos utiliza-la
#tambem iremos utilizar a funcao norma (para verificar o comprimento)

#Entrada: conjunto de vetores
#Etapa 1: verificar se entrada saõ vetores do mesmo tamanho, senao retorna erro
#Etapa 2; verifica se comprimento dos vetores é 1, senao retorna false
#Etapa 2: faz pares de vetores
#Etapa 3: verifica se os pares sao ortogonais, senao retorna false
#Saida; true(é base ortonormal) ou false (nao é base ortonormal)
#Autora: gabriella radke
#-------------------------------------------------------------------------------------
function base_ortonormal(set_of_vectors)
    #checando se o conjunto esta em forma de array
    if isa(set_of_vectors,Array) == true
        #pega tamanho do primeiro vetor, iremos comparar com o resto
        numberOfLines, numberOfColumns = size(set_of_vectors[1])
        #alguma das dimensoes precisa ser 1, senao nao é vetor...
        if numberOfLines != 1 && numberOfColumns != 1
            #se nenhuma for 1, retorna erro
            throw(ArgumentError("elements of set_of_vectors should be vectors (dimensions 1xn or nx1)"))
        else
            #para cada vetor, ele deve ter o mesmo tamanho do primeiro e comprimento 1
            for vector in set_of_vectors
                if numberOfLines == size(vector)[1] && numberOfColumns == size(vector)[2] && 
                    norm(vector) <= 1.001 && norm(vector) >= 0.999
                    continue
                else
                    return false
                # cada par de vetor devera ser ortogonal entre si
                for vector2 in set_of_vectors
                    #if dois_vetores_ortogonais(vector,vector2) == true
                    if dot(vector,vector2) == 0   
                        continue
                    else
                        return false     
                    end
                end
                end
            end
        end
    else 
        #retorna erro se entrada nao for array
        throw(ArgumentError(" 'set_of_vectors' parameter should be an array"))
    end
    #se chegamos ate aqui, é pq a base é ortonormal
    return true
end

#teste para base ortonormal
set_of_vectors = [[1/sqrt(2) 0 -1/sqrt(2)], [1/2 sqrt(2)/2 1/2], [1/2 -sqrt(2)/2 1/2]]

#teste
base_ortonormal([[1/sqrt(2) 0 -1/sqrt(2)], [1/2 sqrt(2)/2 1/2], [1/2 -sqrt(2)/2 1/2]])


#******************************************************************************************
# PROJECAO RETA Calcula a projeção de um vetor em uma reta
# Entradas: um vetor "v" e uma reta "R" (matriz)
# Saídas: float
#------------------------------------------------------------------------------------------
function projecao_vetor_reta(v, R)
    x1, x2 = R; # x1 receberá o vetor do início do segmento da reta "R", e x2 receberá o vetor do final
    return produto_interno(v, x2 - x1) / (norma(x2 - x1) ^ 2) * (x2 - x1);
end


#******************************************************************************************
# PROJECAO PLANO Calcula a projeção de um vetor em um plano
# Entradas: 2 vetores, 1 é o vetor a ser projetado e o segundo é ortogonal ao plano
# Saídas: float
#------------------------------------------------------------------------------------------
function projecao_vetor_plano(v, u)
    return v - produto_interno(v, u) / (norma(u) ^ 2) * u;
end


#******************************************************************************************
# DECOMPOSICAO QR Decomposição QR de uma matriz A.
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


#******************************************************************************************
# DECOMPOSICAO QR Fatoração QR baseado em Projeções
#------------------------------------------------------------------------------------------
function fatorizacao_qr(A)
   m,n = size(A);
   Q = zeros(n,n);
   R = zeros(n,n);
   Q[:,1]=(1/norm(A[:,1]))*A[:,1];
   for j = 2:n
       #Obter Q por Normalizacao do vetor D
       P = zeros(n);
       for k = 1:(j-1)
           P  = P + (A[:,j]'*Q[:,k])*Q[:,k];
       end
       D = A[:,j] - P;
       Q[:,j]=(1/norm(D))*D;
   end
   R = Q'A;
   return Q, R
end



#******************************************************************************************
#TRIANGULAR - verificar se uma dada matriz é triangular
#Notr que ser triangular (inferior ou superior) é uma propridade de matrizes quadradas
#portanto, devemos primeiramente verificar se o input é uma matriz quadrada
#para então olhar suas entradas e caracteriza-la como triangular ou não
#uma matriz é triangular superior se todas as entradas abaixo da diagonal são 0
#é chamada triangular inferior se todas as entradas acima da diagonal são 0.

#Entrada: matriz, se teste é para superior/inferior (default é superior, para inferior usar false em upper)
#Etapa 1: verifica se é matriz (senão retorna erro); 
#Etapa 2: verifica se matriz é quadrada (senão retorna erro);
#Etapa 3: olha elementos abaixo ou acima da diagonal para determinar se é triangular ou não
#Saída: True (é triangular) ou False ( não é triangular)
#Autora: gabriella radke
#-------------------------------------------------------------------------------------

function eh_triangular(matrix, upper = true)

    #checamos se a entrada de fato é uma matriz
    if isa(matrix, Array{Float64,2}) == true  || isa(matrix, Array{Int64,2}) == true
        #checamos se é quadrada 
        numberOfLines, numberOfColumns = size(matrix)
        if numberOfLines == numberOfColumns 
            #se queremos testar se é triangular superior
            if upper == true 
                #checamos se os elementos abaixo da diagonal sao 0 (a_ij i>j)
                for line in 2:numberOfLines
                    for column in 1:line-1 
                        if matrix[line,column] == 0 
                            continue
                        else 
                            #se algum elemento abaixo da diagonal não for 0, nao é triangular superior
                            return false 
                        end
                    end
                end
            #se queremos testar se é triangular inferior
            elseif upper == false 
                #checamos se os elementos acima da diagonal sao 0 (a_ij i<j)
                for line in 1:numberOfLines-1 
                    for column in line+1:numberOfColumns 
                        if matrix[line,column] == 0 
                            continue
                        else 
                            #se algum elemento acima da diagonal não for 0, nao é triangular inferior
                            return false 
                        end
                    end
                end
            else 
                #se o parametro upper nao for true nem false retornamos erro
                throw(ArgumentError(" 'upper' parameter should be true or false"))
            end
        else 
            #retornarmos erro se matriz não for quadrada
            throw(ArgumentError(" 'matrix' parameter should be a SQUARE matrix"))
        end
    
    else 
        #retornamos erro se input nao for array
        throw(ArgumentError("'matrix' parameter should be a matrix")) 
    end    
    #se chegamos até aqui, é pq a bendita é triangular
    return true       
end



#teste
eh_triangular(upperTriangularMatrix)


#******************************************************************************************
# SUBTITUICAO Substituição reversa em uma matriz triangular superior quadrada não-singular.
#
# Entrada: uma matriz quadrada triangular não singular A e um vetor linha b
# Saída: o vetor x, solução de Ax=b.
# Autor: Gastão.
#------------------------------------------------------------------------------------------
function retrosubstituicao_triangular_superior_quadrada(A,b)
    n=length(b)
    x=zeros(n)
    D=0
     
    for i=n:-1:1                             # Indo da última incógnita para a primeira.

        D=A[i,:]'*x                          # Substituímos os valores anteriores.
        
        x[i]=(b[i]-D)/A[i,i]                 # Valor da i-ésima incógnita.
    end
    
    return x 
end


#******************************************************************************************
#SUBSTITUICAO - solucionar sistemas de equaçoes lineares Ax=b 
#a matriz devera ser quadrada e triangular

#Entrada: matriz quadrada A e vetor b
#Etapa 1: verificar se matriz é quadrada nxn vetor tem dimensao correta nx1
#Etapa 2: verificar se matriz é triangular superior ou inferior
#Etapa 3: substituiçao
#Saida: vetor x - solucao do sistema linear
#Autora: gabriella radke
#-------------------------------------------------------------------------------------
function substituicao(matrix,vector) 
    #checa se input é matriz
    if isa(matrix, Array{Float64,2}) == true  || isa(matrix, Array{Int64,2}) == true
        #Se não for quadrada nxn ou vetor nao tiver dimensao correta nx1 retorna erro
        numberOfLines,numberOfColumns = size(matrix)
        if numberOfLines != numberOfColumns || numberOfLines != size(vector)[1] || 1 != size(vector)[2]
            throw(ArgumentError("'matrix' parameter should be a nxn square matrix and 
                    'vector' should be a nx1 array"))
        end
        #caso matriz seja triangular inferior
        if eh_triangular(matrix,false) == true  
            #susbtituicao!
            x = zeros(numberOfLines)
            x[1] = vector[1]/matrix[1,1]
            for i in 2:numberOfLines 
            x[i] = (vector[i]- sum(matrix[i,j]*x[j] for j in 1:i-1))/matrix[i,i]
                        end 
            return x
        #Caso matriz for triangular superior
        elseif eh_triangular(matrix,true) == true
                #substituicao!
                x = zeros(numberOfLines) 
                x[numberOfLines] = vector[numberOfLines]/matrix[numberOfLines,numberOfColumns]
                for i in numberOfLines-1:-1:1
                    x[i] = (vector[i]- sum(matrix[i,j]*x[j] for j in i+1:n))/matrix[i,i]
                                end            
                return x
         #retorna erro se input não é array
        throw(ArgumentError("'matrix' parameter should be a matrix"))
        end
    end
end 

#teste
b = randn(3,1)
A = [1.0      0.0      0.0;
 2.12374  1.0      0.0;
 1.20021  1.04731  1.0]

print(A*substituicao(A,b), b)

#******************************************************************************************
# SUBSTITUICAO Obtenção do vetor X de soluções do Sistema Linear AX = B
# por Substituicao Reversa usando a Matriz Triangular Superior
#-------------------------------------------------------------------------------------
function substituicaoReversa_TS(A,B)
   m,n = size(A);
   X = zeros(n);
   for i = n:-1:1
       S = 0;
       for j = i+1:n
           S = S + A[i,j]*X[j];
       end
       X[i] = (B[i] - S)/A[i,i];
   end
   return X;
end

#******************************************************************************************
# SUBSTITUICAO Obtenção do vetor X de soluções do Sistema Linear AX = B
# por Substituicao Direta usando a Matriz Triangular Inferior
#-------------------------------------------------------------------------------------
function substituicaoDireta_TI(A,B)
   m,n = size(A);
   X = zeros(n);
   for i = 1:n
       S = 0;
       for j = 1:i-1
           S = S + A[i,j]*X[j];
       end
       X[i] = (B[i] - S)/A[i,i];
   end
   return X;
end


#******************************************************************************************
# ELIMINACAO GAUSS Realiza a eliminação gaussiana, usa uma sequência de operações de linha elementares para 
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


#******************************************************************************************
# ELIMICACAO GAUSS Eliminação gaussiana com pivoteamento parcial de uma matriz quadrada não singular com pivôs não-nulos.
#
# Entrada: uma matriz com os pontos dados.
# Saída: o vetor direção "v" da reta que melhor se aproxima dos pontos dados.
# Autor: Gastão.
#------------------------------------------------------------------------------------------
function resolver_escalonamento_com_pivoteamento(A,b)
    A,b=escalonamento_com_pivoteamento(A,b)                        # Escalona o sistema.
    x=retrosubstituicao_triangular_superior_quadrada(A,b)          # Dá a solução do sistema pelo escalonamento.
   return x
end


#******************************************************************************************
# ELIMINICAO GAUSS Processo de Transformação das matrizes A e B 
# para obter uma Matriz Triangular Superior
#------------------------------------------------------------------------------------------
function eliminacaoGaussiana_TS(A,B)
   m,n = size(A);
   for j = 1:n-1
       for i = j+1:n
           E = -1*(A[i,j]/A[j,j]);
           A[i,:] = A[i,:] + E*A[j,:];
           B[i]   = B[i]   + E*B[j];
       end
   end
   return A,B
end

#******************************************************************************************
# ELIMINACAO GAUSS Processo de Transformação das matrizes A e B 
# para obter uma Matriz Triangular Inferior
#------------------------------------------------------------------------------------------
function eliminacaoGaussiana_TI(A,B)
   m,n = size(A);
   for j = n:-1:2
       for i = 1:j-1
           E = -1*(A[i,j]/A[j,j]);
           A[i,:] = A[i,:] + E*A[j,:];
           B[i]   = B[i]   + E*B[j];
       end
   end
   return A,B
end


#******************************************************************************************
# GAUSS JORDAN Processo de Transformação das matrizes A e B 
# para obter uma Matriz Diagonal
#------------------------------------------------------------------------------------------
function GaussJordan(A,B)
   m,n = size(A);
   X = zeros(n);
   for j = 1:n
       for i = 1:n
           if i!=j
              E = -1*(A[i,j]/A[j,j]);
              A[i,:] = A[i,:] + E*A[j,:];
              B[i]   = B[i]   + E*B[j];
           end
       end
       X[j] = B[j]/A[j,j];
   end
   return X;
end

#******************************************************************************************
#  DECOMPOSICAO L1U Fatoração L1U baseado em Aproximaçôes de Matrizes de Posto 1
# Diagonal da Matriz Triangular Inferior L con 1s
#------------------------------------------------------------------------------------------
function fatoracao_L1U(A)
  m,n = size(A);
  p = min(m,n);
  L = zeros(m,p);
  U = zeros(p,n);
  T = copy(A);
  for k = 1:p
      L[:,k] = (1/T[k,k])*T[:,k];
      U[k,:] = T[k,:];
      Z = L[:,k]*U[k,:]';
      T = T - Z;
  end
  return L,U
end

#******************************************************************************************
# DECOMPOSICAO LU Fatoração LU1 baseado em Aproximaçôes de Matrizes de Posto 1
# Diagonal da Matriz Triangular Superior U con 1s
#------------------------------------------------------------------------------------------
function Factorization_LU1(A)
  m,n = size(A);
  p = min(m,n);
  L = zeros(m,p);
  U = zeros(p,n);
  T = copy(A);
  for k = 1:p
      L[:,k] = T[:,k];
      U[k,:] = (1/T[k,k])*T[k,:];
      Z = L[:,k]*U[k,:]';
      T = T - Z;
  end
  return L,U
end
#----------------------------------------------------------------------
function LU(A)
  return fatoracao_L1U(A);
end
function LU(A,Sw)
  if Sw
     return fatoracao_L1U(A);
  else
     return Factorization_LU1(A);
  end 
end



#******************************************************************************************
# POSTO Calcula o posto de uma Matriz, que corresponte ao número de linhas ou colunas 
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

#******************************************************************************************
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


#******************************************************************************************
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


#******************************************************************************************
# ESCALONAMENTO Escalonamento com pivoteamento parcial de uma matriz quadrada não singular com pivôs não-nulos.
#
# Entrada: uma matriz quadrada não singular e um vetor linha b.
# Saída: uma matriz quadrada não singular escalonada e um vetor linha b que geram
# um sistema equivalente ao dos dados de entrada.
# Autor: Gastão.
#------------------------------------------------------------------------------------------
function escalonamento_com_pivoteamento(A,b)
    A_amp=[A b']                                   # Matriz ampliada de A
    n,m=size(A)
  
    for pivo=1:(n-1)                                
        A_amp=pivoteamento_parcial(A_amp,pivo)     # Encontra um "novo" pivô.
        A_amp=zerar_abaixo(A_amp,pivo)             # Zera os elementos abaixo do pivô, 
    end
    
    A=A_amp[:,1:m]                                 # Matrizes resultantes do escalonamento.
    b=A_amp[:,m+1]
    return A,b
end
#------------------------------------------------------------------------------------------
function pivoteamento_parcial(A,pivo)
    Ap=copy(A)
    n_pivo=novo_pivo(Ap[:,pivo],pivo)              # "Novo" pivô (pode ser o original).
    
    linha=Ap[pivo,:]                               # Troca as linhas do pivô original com a do novo pivô.
    Ap[pivo,:]=Ap[n_pivo,:]
    Ap[n_pivo,:]=linha
    
    return Ap
end
#------------------------------------------------------------------------------------------
function novo_pivo(coluna,pivo)
    
    maximo=abs(coluna[pivo])             # Pivô original
    n_pivo=pivo
        
    for i=pivo:length(coluna)            # Procura o maior valor entre as células abaixo da célula
        if (abs(coluna[i]))>maximo       # do pivô e transforma indica essa célula como a do novo pivô.
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
        A[i,:]+=-(A[i,pivo]/A[pivo,pivo])A[pivo,:]     # Modifica todas as linhas abaixo da linha do pivô,
    end                                                # zerando os elementos abaixo deste.
    
    return A
end

#******************************************************************************************
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

#******************************************************************************************
# JACOBI Metodo de Jacobi
# Realiza Aproximações Sucesivas para obter a Solução
# Usa o vetor auxiliar Y para guardar a solucao da iteracao anterior
#------------------------------------------------------------------------------------------
function jacobi(A,B,E)
   m,n = size(A);
   Y = rand(n);
   X = zeros(n);
   t = 0;
   while NormVector(X-Y)>E
      t = t + 1;
      Y = copy(X)
      for i = 1:n
          S = 0;
          for j = 1:n
              if i!=j
                 S = S + A[i,j]*Y[j];
              end
          end
          X[i] = (B[i]- S)/A[i,i];
      end
   end
   return X
end

#******************************************************************************************
# GAUSS SEIDEL Metodo de Gauss-Seidel
# Realiza Aproximações Sucesivas para obter a Solução
# Usa o vetor auxiliar Y para guardar a solucao da iteracao anterior
#------------------------------------------------------------------------------------------
function gauss_Seidel(A,B,E)
   m,n = size(A);
   Y = rand(n);
   X = zeros(n);
   t = 0;
   while NormVector(X-Y)>E
      t = t + 1;
      Y = copy(X)
      for i = 1:n
          S = 0;
          for j = 1:n
              if i!=j
                 S = S + A[i,j]*X[j];
              end
          end
          X[i] = (B[i]- S)/A[i,i];
      end
   end
   return X
end



#******************************************************************************************
# F U N Ç Õ E S   D E   V I S U A L I Z A Ç Ã O
#******************************************************************************************
#******************************************************************************************
# Formatea um Número Real em "PI" Posições Inteiras e 
# "PD" Posições Decimais
#------------------------------------------------------------------------------------------
function formatarReal(X,PI,PD)
  X = round(X; digits=PD);
  S = string(X);
  L = length(S);
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SwNegative = (S[1]=='-');
  if SwNegative==true
     if X==0
        S = "0.0";
        L = 3;
        SwNegative = false;
     else
        S = S[2:L];
        L = L - 1;
     end
  end
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  n = 0;
  Q = 0;
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  E = findfirst('e', S);
  if E==nothing
     E = 0;
  else
     Q = findfirst('-', S);
     if Q==nothing
        Q = 0;
     end
     n = parse(Int,S[(Q==0 ? E+1 : Q+1):L]);  
    #n = (Q==0 ? parse(Int,S[E+1:L]) : parse(Int,S[Q+1:L]));
     S = S[1:E-1];
     L = E-1;
  end
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  P = findfirst('.', S);
  if P==nothing
     S = string(S,".0");
     L = L + 2;
     P = L - 1;
  end
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  NI = S[1:P-1];
  LI = P - 1;
  ND = S[P+1:L];
  LD = L - P;
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if E==0
  else
     if Q==0
        ND = string(ND, "0"^(n<LD ? 0 : n-LD+1));
        LD = length(ND);
     else
        NI = string("0"^(n<LI ? 0 : n-LI+1), NI);
        LI = length(NI);
     end
     WW = NI;
     NI = (Q==0 ? string(NI,ND[1:n]) : NI[1:LI-n]);
     ND = (Q==0 ? ND[n+1:LD]         : string(WW[LI-n+1:LI],ND));
     LI = length(NI);
     LD = length(ND);
  end
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  NI = (SwNegative==true ? string("-",NI) : NI );
  LI = (SwNegative==true ? LI+1 : LI );
  NI = (PI<=LI ? NI       : string(" "^(PI-LI),NI) );
  ND = (PD<=LD ? ND[1:PD] : string(ND,"0"^(PD-LD)) );
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  return string(NI,".",ND);
end

#******************************************************************************************
# Imprime Horizontalmente o vetor U considerando 
# "PI" Posições Inteiras e "PD" Posições Decimais
#------------------------------------------------------------------------------------------
function imprimeVetorH(U,PI,PD)
   n = length(U)
   for k = 1:n
       print(formatarReal(U[k],PI,PD));
   end
end

#******************************************************************************************
# Imprime Verticalmente o vetor U considerando 
# "PI" Posições Inteiras e "PD" Posições Decimais
#------------------------------------------------------------------------------------------
function imprimeVetorV(U,PI,PD)
   n = length(U)
   for k = 1:n
       println(formatarReal(U[k],PI,PD));
   end
end

#******************************************************************************************
# Imprime os elementos da matriz A considerando 
# "PI" Posições Inteiras e "PD" Posições Decimais
#------------------------------------------------------------------------------------------
function imprimeMatriz(A,PI,PD)
   m,n = size(A)
   for i = 1:m
       for j = 1:n
           print(formatarReal(A[i,j],PI,PD));
       end
       println("");
   end
end

#******************************************************************************************
# Imprime os elementos da matriz A e do vetor B 
# considerando "PI" Posições Inteiras e "PD" Posições Decimais
#------------------------------------------------------------------------------------------
function imprimeSistemaLinear(A,B,PI,PD)
   m,n = size(A)
   for i = 1:m
       for j = 1:n
           print(formatarReal(A[i,j],PI,PD));
       end
       print(" "^PI,"│");
       print(formatarReal(B[i],PI,PD));
       println("");
   end
end

#******************************************************************************************
# Imprime a Matriz em um Arquivo TXT
#------------------------------------------------------------------------------------------
function matrizParaTexto(A,PI,PD,FILE)
   m,n = size(A);
   F = open(FILE,"w");
   for i = 1:m
       for j = 1:n
           write(F,formatarReal(A[i,j],PI,PD));
       end
       write(F,"\n");
   end
   close(F);
end
#------------------------------------------------------------------------------------------
#******************************************************************************************

#******************************************************************************************
# Gera Codigo CSS para formato da Pagina HTML
#------------------------------------------------------------------------------------------
function pegarCSS()
   ForeColorCell     = "#000000";  #Text Color of Cell
   BackColorCell     = "#FFDB4D";  #Background Color of Cell
   ForeColorIndex    = "#FFFFFF";  #Text Color of Index
   BackColorIndex    = "#505050";  #Background Color of Index
   ForeColorDiagonal = "#000000";  #Text Color of Diagonal
   BackColorDiagonal = "#99E6FF";  #Background Color of Diagonal
   TitleColor        = "#0000FF";  #Title Color
   TextColor         = "#00FF00";  #Text Color
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   S = "<STYLE>\n";
   S = string(S,".TABLE     {  padding: 10px; border-spacing:1px 1px; }\n" );
   S = string(S,".INDEX     {  color: " , ForeColorIndex    , "; background-color: " , BackColorIndex    , "; text-align: center;  padding: 10px;  font-family: Calibri, Tahoma, Verdana; }\n" );
   S = string(S,".DIAGONAL  {  color: " , ForeColorDiagonal , "; background-color: " , BackColorDiagonal , "; text-align: center;  padding: 15px;  font-family: 'Courier New', Tahoma, Verdana; font-weight: bold; }\n" );
   S = string(S,".CELL      {  color: " , ForeColorCell     , "; background-color: " , BackColorCell     , "; text-align: center;  padding: 15px;  font-family: 'Courier New', Tahoma, Verdana; }\n" );
   S = string(S,".TITLE     {  color: " , TitleColor        , "; text-align: left;  padding: 15px;  font-family: Tahoma, Verdana, Calibri;  font-size: 28px; }\n" );
   S = string(S,".TEXT      {  color: " , TextColor         , "; text-align: left;  padding: 15px;  font-family: Tahoma, Verdana, Calibri;  font-size: 28px; }\n</STYLE>\n" );
   S = string(S,"</STYLE>\n" );
   #-------------------------------------------------------------------------------
   return string(S);
end

#******************************************************************************************
# Gera o Codigo HTML associado para a Matriz 
#------------------------------------------------------------------------------------------
function escreveMatriz(F,A,PD,ShowIndex,ShowDiagonal)
   m,n = size(A);
   #-------------------------------------------------------------------------------
   write(F,"<TABLE class=TABLE>\n");
   if ShowIndex==true
      write(F,"<TR><TD> </TD>");
      for j = 1:n
          write(F,"<TD CLASS=INDEX>",string(j),"</TD>" );
      end
      write(F,"</TR>\n");
   end
   for i = 1:m
       write(F,"<TR ALIGN=CENTER>");
       if ShowIndex==true
          write(F,"<TD CLASS=INDEX>",string(i),"</TD>" );
       end
       for j = 1:n
           write(F, (ShowDiagonal==true && i==j ? "<TD CLASS=DIAGONAL>" : "<TD CLASS=CELL>"), formatarReal(A[i,j],0,PD), "</TD>" );
       end
       write(F,"</TR>\n");
   end
   write(F,"</TABLE><br><br>\n</BODY>\n</HTML>\n");
end

#******************************************************************************************
# Resultados da Matriz em Arquivo HTML
#------------------------------------------------------------------------------------------
function MatrixToHTML(A,PD,FILE,ShowIndex,ShowDiagonal,Title)
   #-------------------------------------------------------------------------------
   F = open(FILE,"w");
   #-------------------------------------------------------------------------------
   write(F,"<HTML>\n<HEAD>\n");
   write(F,pegarCSS());
   write(F,"</HEAD>\n<BODY>\n");
   #-------------------------------------------------------------------------------
   write(F,"<H1>",Title,"</H1>\n","<HR>\n");
   #-------------------------------------------------------------------------------
   escreveMatriz(F,A,PD,ShowIndex,ShowDiagonal);
   #-------------------------------------------------------------------------------
   close(F);
end

#******************************************************************************************
# Simulação da Fatoração QR
#------------------------------------------------------------------------------------------
function simulacao_QR(A,PD,FILE)
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   F = open(FILE,"w");
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write(F,"<HTML>\n<HEAD>\n");
   write(F,pegarCSS());
   write(F,"</HEAD>\n<BODY>\n");
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write(F,"<H1>Simula&ccedil;&atilde;o do M&eacute;todo QR</H1>\n" );
   write(F,"<HR>\n");
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write(F,"<SPAN CLASS=TITLE>Matriz A</SPAN>\n");
   escreveMatriz(F,A,PD,true,false);
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   m,n = size(A);
   Q = zeros(n,n);
   R = zeros(n,n);
   Q[:,1]=(1/norm(A[:,1]))*A[:,1];
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write(F,"<SPAN CLASS=TITLE>Passo 1: Inicio</SPAN>\n<HR>\n");
   write(F,"<SPAN CLASS=TITLE>Matriz Q</SPAN>\n");
   escreveMatriz(F,Q,PD,true,false);
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   for j = 2:n
       P = zeros(n)
       for k = 1:(j-1)
           P  = P + (A[:,j]'*Q[:,k])*Q[:,k]
       end         
       D = A[:,j] - P
       Q[:,j]=(1/norm(D))*D
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       write(F,"<SPAN CLASS=TITLE>Passo ",string(j),"</SPAN>\n<HR>\n");
       write(F,"<SPAN CLASS=TITLE>Matriz Q</SPAN>\n");
       escreveMatriz(F,Q,PD,true,false);
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   end
   R = Q'A 
   write(F,"</BODY>\n</HTML>\n");
   close(F);
end

#******************************************************************************************
# Simulação da Fatoração L1U por Aproximação de Matriz de Posto 1
#------------------------------------------------------------------------------------------
function simulacao_L1U(A,PD,FILE)
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   F = open(FILE,"w");
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write(F,"<HTML>\n<HEAD>\n");
   write(F,pegarCSS());
   write(F,"</HEAD>\n<BODY>\n");
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write(F,"<H1>Simula&ccedil;&atilde;o do M&eacute;todo L1U por Aproxima&ccedil;&atilde;o de Matriz de Posto 1</H1>\n" );
   write(F,"<HR>\n");
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write(F,"<SPAN CLASS=TITLE>Matriz A</SPAN>\n");
   escreveMatriz(F,A,PD,true,true);
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   m,n = size(A);
   p = min(m,n);
   L = zeros(m,p);
   U = zeros(p,n);
   T = copy(A);
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write(F,"<SPAN CLASS=TITLE>Passo 0: Inicio</SPAN>\n<HR>\n");
   write(F,"<SPAN CLASS=TITLE>Matriz L</SPAN>\n");
   escreveMatriz(F,L,PD,true,true);
   write(F,"<SPAN CLASS=TITLE>Matriz U</SPAN>\n");
   escreveMatriz(F,U,PD,true,true);
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   for k = 1:p
       L[:,k] = (1/T[k,k])*T[:,k];
       U[k,:] = T[k,:];
       Z = L[:,k]*U[k,:]';
       T = T - Z;
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       write(F,"<SPAN CLASS=TITLE>Passo ",string(k),"</SPAN>\n<HR>\n");
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       write(F,"<SPAN CLASS=TITLE>Matriz de Posto 1</SPAN>\n");
       escreveMatriz(F,Z,PD,true,true);
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       write(F,"<SPAN CLASS=TITLE>Matriz Resultante</SPAN>\n");
       escreveMatriz(F,T,PD,true,true);
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       write(F,"<SPAN CLASS=TITLE>Matriz L</SPAN>\n");
       escreveMatriz(F,L,PD,true,true);
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       write(F,"<SPAN CLASS=TITLE>Matriz U</SPAN>\n");
       escreveMatriz(F,U,PD,true,true);
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   end
   write(F,"</BODY>\n</HTML>\n");
   close(F);
end

#******************************************************************************************
# Simulação da Fatoração LU1 por Aproximação de Matriz de Posto 1
#------------------------------------------------------------------------------------------
function simulacao_LU1(A,PD,FILE)
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   F = open(FILE,"w");
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write(F,"<HTML>\n<HEAD>\n");
   write(F,pegarCSS());
   write(F,"</HEAD>\n<BODY>\n");
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write(F,"<H1>Simula&ccedil;&atilde;o do M&eacute;todo LU1 por Aproxima&ccedil;&atilde;o de Matriz de Posto 1</H1>\n" );
   write(F,"<HR>\n");
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write(F,"<SPAN CLASS=TITLE>Matriz A</SPAN>\n");
   escreveMatriz(F,A,PD,true,true);
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   m,n = size(A);
   p = min(m,n);
   L = zeros(m,p);
   U = zeros(p,n);
   T = copy(A);
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write(F,"<SPAN CLASS=TITLE>Passo 0: Inicio</SPAN>\n<HR>\n");
   write(F,"<SPAN CLASS=TITLE>Matriz L</SPAN>\n");
   escreveMatriz(F,L,PD,true,true);
   write(F,"<SPAN CLASS=TITLE>Matriz U</SPAN>\n");
   escreveMatriz(F,U,PD,true,true);
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   for k = 1:p
       L[:,k] = T[:,k];
       U[k,:] = (1/T[k,k])*T[k,:];
       Z = L[:,k]*U[k,:]';
       T = T - Z;
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       write(F,"<SPAN CLASS=TITLE>Passo ",string(k),"</SPAN>\n<HR>\n");
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       write(F,"<SPAN CLASS=TITLE>Matriz de Posto 1</SPAN>\n");
       escreveMatriz(F,Z,PD,true,true);
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       write(F,"<SPAN CLASS=TITLE>Matriz Resultante</SPAN>\n");
       escreveMatriz(F,T,PD,true,true);
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       write(F,"<SPAN CLASS=TITLE>Matriz L</SPAN>\n");
       escreveMatriz(F,L,PD,true,true);
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       write(F,"<SPAN CLASS=TITLE>Matriz U</SPAN>\n");
       escreveMatriz(F,U,PD,true,true);
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   end
   write(F,"</BODY>\n</HTML>\n");
   close(F);
end
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------