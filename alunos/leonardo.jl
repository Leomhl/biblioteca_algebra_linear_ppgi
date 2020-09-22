#------------------------------------------------------------------------------------------
# Funções auxiliares 
#------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
# Retorna o valor máximo de um vetor
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
# Funções vetoriais
#------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
# Produto interno
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
# Norma de um vetor (comprimento)
# Entradas: 1 vetores
# Saídas: 1 valor numérico
#------------------------------------------------------------------------------------------
function norma(a) 
    return sqrt(produto_interno(a, a))
end 


#------------------------------------------------------------------------------------------
# Testa a ortogonalidade entre dois vetores
# Entradas: 2 vetores
# Saídas: booleano (true ou false)
#------------------------------------------------------------------------------------------
function dois_vetores_ortogonais(a, b, tol=10e-5)
    return produto_interno(a, b) <= tol
end


#------------------------------------------------------------------------------------------
# Testa a ortogonalidade entre um vetor e o plano
# Entradas: 1 vetor, 1 matriz que constrói o plano
# Saídas: booleano (true ou false)
#------------------------------------------------------------------------------------------
function vetor_ortogonal_plano(v, P, tol=10e-5)
    return norma(multiplicar_matrizes(P, v)) <= tol
end

