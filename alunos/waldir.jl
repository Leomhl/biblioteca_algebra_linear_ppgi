# funções auxiliares

function max_val(v)
    m = -Inf
    for i in v
        if(i > m) m = i end
    end
    return m
end

function min_val(v)
    m = Inf
    for i in v
        if(i < m) m = i end
    end
    return m
end

function ensure_matrix(v) #se v for um vetor, transforma v em uma matriz de 1 coluna. Se v tiver mais de 2 dimensões, causa um erro
    @assert 1 <= ndims(v) <= 2 
    return v[:,:]
end

function ensure_vector(v) #se v for um matriz de 1 coluna/linha, transforma v em um vetor
    @assert ndims(v) == 1 || (max_val(size(v)) == length(v))
    return v[:]
end

function mat_mul(A, B)
    A = ensure_matrix(A)
    B = ensure_matrix(B)

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

function mat_sum(A, B)
    A = ensure_matrix(A)
    B = ensure_matrix(B)

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


# IX - Se dois vetores são ortogonais


# dois_vetores_ortogonais(a, b): recebe 2 vetores "u" e "v" e retorna se eles são ortogonais. 
function dois_vetores_ortogonais(a, b, tol=10e-5)
    return dot(a, b) <= tol
end

# X - Se uma vetor é ortogonal a um plano
# v é um vetor, e P é u
function vetor_ortogonal_plano(v, P, tol=10e-5)
    return norm(mat_mul(P, v)) <= tol
end

# I - Se duas retas são ortogonais

# retas_ortogonais(a, b): recebe 2 retas "a" e "b" e retorna se elas são ortogonais. Cada reta é representada por 2 vetores.
# "a" e "b" são tuplas com 2 vetores cada uma, que representam o início e o final de segmentos de retas
function retas_ortogonais(a, b)
    x1, x2 = a; # x1 receberá o vetor do início do segmento da reta "a", e x2 receberá o vetor do final
    y1, y2 = b; # y1 receberá o vetor do início do segmento da reta "b", e y2 receberá o vetor do final
    return dois_vetores_ortogonais(x2 - x1, y2 - y1)
end

# IV - Normas

# norma_do_maximo(v): retorna a norma do máximo de v
function norma_do_maximo(a) 
    return max_val(abs.(a))
end 


# II - Projeção de um vetor em uma reta

# projecao_vetor_reta(v, R): recebe um vetor "v" e uma reta "R" e retorna a projeção de v em R
function projecao_vetor_reta(v, R)
    x1, x2 = R; # x1 receberá o vetor do início do segmento da reta "R", e x2 receberá o vetor do final
    return dot(v, x2 - x1) / (norm(x2 - x1) ^ 2) * (x2 - x1);
end

# III - Projeção de um vetor em um plano: recebe um vetor "v" e um vetor "u" ortogonal ao plano P
function projecao_vetor_plano(v, u)
    return v - dot(v, u) / (norm(u) ^ 2) * u;
end

# V - Produto de matrizes (ou recursiva)

# prod_mat_rec(A, B, C, D, ...): recebe um número qualquer de matrizes e retorna o produto entre estas, calculado de forma recursiva
function prod_mat_rec(A, B...) # B... significa um número variável de argumentos
    if(length(B) == 0) # se não houver nenhum argumento adicional, retorna a própria matriz A
        return A
    else # se houver algum argumento adicional, retorna a matriz A multiplicada pelo produto das matrizes restantes
        return mat_mul(A, prod_mat_rec(B...))
    end
end

# VI - Soma de matrizes (ou recursiva)

# soma_mat_rec(A, B, C, D, ...): recebe um número qualquer de matrizes e retorna a soma entre estas, calculada de forma recursiva
function soma_mat_rec(A, B...) # B... significa um número variável de argumentos
    if(length(B) == 0) # se não houver nenhum argumento adicional, retorna a própria matriz A
        return A
    else # se houver algum argumento adicional, retorna a matriz A multiplicada pelo produto das matrizes restantes
        return mat_sum(A, soma_mat_rec(B...))
    end
end

#testes pra ver se as funções estão corretas
function testes()
    # using LinearAlgebra
    tol = 10e-5;

    @assert norma_do_maximo([1,-1]) == 1
    @assert abs(norm([1,-1]) - sqrt(2)) < tol

    r = ( [ 1, 1, 1 ], [ 2, 1, -3 ] );
    s = ( [ 0, 1, 0 ], [-1, 2,  0 ] );
    @assert !retas_ortogonais(r,s)

    v = [ 1, 2, 3 ];
    R = [[ 5, 6, 2 ], [ 0, 0, 0 ]];
    proj = projecao_vetor_reta(v, R);
    esperado = [1.76923077, 2.12307692, 0.70769231];
    @assert norm(proj - esperado) < tol;

    v = [ 2, 5, 8 ];
    u = [ 1, 1, 7];

    proj = projecao_vetor_plano(v, u);
    esperado = [ 0.76470588,  3.76470588, -0.64705882];
    @assert norm(proj - esperado) < tol;

    @assert norma_do_maximo([-1, -3, 7]) == 7;
    @assert norma_do_maximo([-1, -3, -7]) == 7;

    M = rand(5, 5);
    @assert sum(abs.(prod_mat_rec(M,M,M,M) - M*M*M*M)) < tol;
    @assert sum(abs.(soma_mat_rec(M,M,M,M) - (M+M+M+M))) < tol;

    @assert dot(u, v) == u' * v
end

