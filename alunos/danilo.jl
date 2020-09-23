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
# Teste 1
# Verifica se matrizes aleatórias de dimensões diferentes são transpostas corretamente.
# A verificação compara cada linha da matriz original, com a coluna da matriz transposta.
# Se a transposição aconteceu corretamente, elas devem ser iguais
#--------------------------------------------------------------------------------------

function verificacao_transposta(A, A_transposta)
    
    if size(A) == (1,)             #verifica se A é um vetor unitário.
        return A == a_transposta   #Se sim, o compara diretamente com a transposta
    else
    
        m,n = size(A)

        for i=1:m
            if A[i,:] == A_transposta[:,i]

            else
                return false
                break
            end
        end
        return true
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
a = [1]
a_transposta = transposta(a)
println(verificacao_transposta(a, a_transposta))