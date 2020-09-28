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

#Verificando uma matriz aleatória 100x100
a = rand(100,100)
a_transposta = transposta(a)
println(verificacao_transposta(a, a_transposta))

#Verificando um vetor 100x1
a = rand(100,1)
a_transposta = transposta(a)
println(verificacao_transposta(a, a_transposta))

#Verificando um vetor 1x100
a = rand(1,100)
a_transposta = transposta(a)
println(verificacao_transposta(a, a_transposta))

#Verificando um vetor nulo (0)
a = [0]
a_transposta = transposta(a)
println(verificacao_transposta(a, a_transposta))

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


#Verificando uma matriz aleatória 100x100
a = rand(100,100)
a_inversa = inversa(a)
println(testa_matriz_inversa(a, a_inversa))

#Verificando um vetor 1x100
a = rand(1,100)
a_inversa = inversa(a)

#Verificando um vetor nulo (0)
a = [0]
a_transposta = transposta(a)
println(verificacao_transposta(a, a_transposta))
