#-----------------------------------------------------------------------------------
#Funções a serem implementadas:
#Triangular - verificar se uma dada matriz é triangular
#Ortogonal - verificar se uma dada matriz é ortogonal
#Ortonormal - verificar se um dado conjunto de vetores forma base ortonormal
#Substituição reversa - para solução de sistemas de equações lineares
#Pseudo inversa - utilizando decomposição svd
#SVD - decomposição (para poder fazer pseudo inversa)
#-----------------------------------------------------------------------------------

using LinearAlgebra

#matrizes para testes
notTriangularMatrix = [-1 1/4 1/4 0;
                        1/4 -1 0 1/4;
                        1/4 0 -1 1/4;
                        0 1/4 1/4 -1]

squareTriU = [-1 1/4 1/4 0;
              0 -1 -1 1/4;
              0 0 1/4 1/4;
              0 0 0 -1]

squareTriL = [-1 0 0 0;
              1/4 -1 0 0;
              1/4 0 -1 0;
              0 1/4 1/4 -1]

rectangularTriU = [ 1 3 1;
                    0 1 7;
                    0 0 4;
                    0 0 0]

rectangularTriL = [2 0 0 0 0; 
                   3 1 0 0 0;
                   4 6 2 0 0]



matriz = lowerTriangularMatrix'

for i in 1:5 
    print(i)
end

#-----------------------------------------------------------------------------------
#Triangular - verificar se uma dada matriz é triangular
#é chamada triangular inferior se todas as entradas acima da diagonal principal são 0.
#é chamada triangular superior se todas as entradas abaixo da diagonal principal são 0

#Entrada: matriz, se teste é para superior/inferior (default é superior, para inferior usar false em upper)
#Etapa 1: transforma entrada em array (caso seja de outro tipo); 
#Etapa 2: olha elementos abaixo ou acima da diagonal principal para determinar se é triangular ou não
#Saída: True (é triangular) ou False ( não é triangular)

function e_triangular(matriz, superior = true)

    #Transformamos entrada em array
    matriz = Array(matriz)
        #pegamos dimensoes da matriz
        numeroDeLinhas, numeroDeColunas = size(matriz)
            #se queremos testar se é triangular superior
            if superior == true 
                #checamos se os elementos abaixo da diagonal sao 0 (a_ij i>j)
                for linha = 2:numeroDeLinhas, coluna = 1:linha-1 
                        #se algum elemento abaixo da diagonal não for 0, nao é triangular superior
                        if matriz[linha,coluna] != 0 
                            return false 
                        end
                end
            #se queremos testar se é triangular inferior
            elseif superior == false 
                #checamos se os elementos acima da diagonal sao 0 (a_ij i<j)
                for linha = 1:numeroDeLinhas-1, coluna = linha+1:numeroDeColunas 
                        #se algum elemento acima da diagonal não for 0, nao é triangular inferior
                        if matriz[linha,coluna] != 0 
                            return false 
                        end
                end
            else 
                #se o parametro upper nao for true nem false retornamos erro
                #so vai retornar esse erro se o parametro superior for uma outra variavel definida 
                #que nao seja true ou false
                throw(ArgumentError("Parâmetro 'superior' deve ser true ou false. Para testar se é
                é triangular inferior, utilizar false."))
            end  
    #se chegamos até aqui, é pq a bendita é triangular
    return true       
end



#teste
e_triangular(notTriangularMatrix, Array(matriz))

#-----------------------------------------------------------------------------------
#Transposta 
#matriz tal que o elemento a_ij é igual ao elemento a_ji da matriz original
#necessária para verificar ortogonalidade

#Entrada: matriz (qualquer tamanho)
#Etapa 1: transforma entrada em array
#Etapa 2: construir matriz transposta vazia com dimensoes corretas
#Etapa 3: popula matriz transposta onde o novo elemento a_ij é igual ao elemento a_ji da matriz original
#Saída: matriz (transposta)

function transposta(matriz)
    
    #trasnforma entrada em array
    matriz = Array(matriz)
    #dimensoes da matriz
    numeroDeLinhas, numeroDeColunas = size(matriz)
    #cria matriz transposta vazia com dimensoes corretas para popularmos depois
    matrizTransposta = zeros(numeroDeColunas, numeroDeLinhas)
    #o elemento a_ij da nova matriz sera o a_ji da original
    for linha = 1:numeroDeColunas, coluna = 1:numeroDeLinhas
        matrizTransposta[linha,coluna] = matriz[coluna,linha]
    end  
    return matrizTransposta
end

#teste
e_triangular(transposta(squareTriL))
#transposta(squareTriU)

#-----------------------------------------------------------------------------------
#Ortogonal - verificar se uma matriz é ortogonal
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

function orthogonal(matriz)
     #checa se input é matriz
    if isa(matriz, Array{Float64,2}) == true  || isa(matriz, Array{Int64,2}) == true
        #checa se input é matriz quadrada
        numeroDeLinhas, numeroDeColunas = size(matriz)
        if numeroDeLinhas == numeroDeColunas
            matrizTransposta = transposta(matriz)
            #se matrix^T * matriz = Id, é ortogonal
            if matrizTransposta*matriz == I
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
orthogonal([1 0; 0 1])

#-----------------------------------------------------------------------------------
#Ortonormal - verificar se um dado conjunto de vetores forma base ortonormal
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

function orthonormal_base(set_of_vectors)
    #checando se o conjunto esta em forma de array
    if isa(set_of_vectors,Array) == true
        #pega tamanho do primeiro vetor, iremos comparar com o resto
        numeroDeLinhas, numeroDeColunas = size(set_of_vectors[1])
        #alguma das dimensoes precisa ser 1, senao nao é vetor...
        if numeroDeLinhas != 1 && numeroDeColunas != 1
            #se nenhuma for 1, retorna erro
            throw(ArgumentError("elements of set_of_vectors should be vectors (dimensions 1xn or nx1)"))
        else
            #para cada vetor, ele deve ter o mesmo tamanho do primeiro e comprimento 1
            for vector in set_of_vectors
                if numeroDeLinhas == size(vector)[1] && numeroDeColunas == size(vector)[2] && 
                    norm(vector) <= 1.001 && norm(vector) >= 0.999
                    continue
                else
                    return false
                #cada par de vetor devera ser ortogonal entre si
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
orthonormal_base([[1/sqrt(2) 0 -1/sqrt(2)], [1/2 sqrt(2)/2 1/2], [1/2 -sqrt(2)/2 1/2]])

#-----------------------------------------------------------------------------------
#Substituição - solucionar sistemas de equaçoes lineares Ax=b 
#a matriz devera ser quadrada e triangular

#Entrada: matriz quadrada A e vetor b
#Etapa 1: verificar se matriz é quadrada nxn vetor tem dimensao correta nx1
#Etapa 2: verificar se matriz é triangular superior ou inferior
#Etapa 3: substituiçao
#Saida: vetor x - solucao do sistema linear

function substitution(matriz,vector) 
    #checa se input é matriz
    if isa(matriz, Array{Float64,2}) == true  || isa(matriz, Array{Int64,2}) == true
        #Se não for quadrada nxn ou vetor nao tiver dimensao correta nx1 retorna erro
        numeroDeLinhas,numeroDeColunas = size(matriz)
        if numeroDeLinhas != numeroDeColunas || numeroDeLinhas != size(vector)[1] || 1 != size(vector)[2]
            throw(ArgumentError("'matrix' parameter should be a nxn square matrix and 
                    'vector' should be a nx1 array"))
        end
        #caso matriz seja triangular inferior
        if is_triangular(matriz,false) == true  
            #susbtituicao!
            x = zeros(numeroDeLinhas)
            x[1] = vector[1]/matriz[1,1]
            for i in 2:numeroDeLinhas 
            x[i] = (vector[i]- sum(matriz[i,j]*x[j] for j in 1:i-1))/matriz[i,i]
                        end 
            return x
        #Caso matriz for triangular superior
        elseif is_triangular(matriz,true) == true
                #substituicao!
                x = zeros(numeroDeLinhas) 
                x[numeroDeLinhas] = vector[numeroDeLinhas]/matriz[numeroDeLinhas,numeroDeColunas]
                for i in numeroDeLinhas-1:-1:1
                    x[i] = (vector[i]- sum(matriz[i,j]*x[j] for j in i+1:n))/matriz[i,i]
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

print(A*substitution(A,b), b)

