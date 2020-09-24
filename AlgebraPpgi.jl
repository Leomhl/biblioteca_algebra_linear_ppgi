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
# Calcula o posto de uma Matriz, que corresponte ao número de linhas ou colunas 
# linearmente independentes da matriz.
# 
# Entradas: 1 vetores
# Saídas: inteiro
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