#-----------------------------------------------------------------------------------
#Funções a serem implementadas:
#Triangular - verificar se uma dada matriz é triangular
#Transposta - calcular matriz transposta
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

upperTriangularMatrix = [-1 1/4 1/4 0;
                        0 -1 -1 1/4;
                        0 0 1/4 1/4;
                        0 0 0 -1]

lowerTriangularMatrix = [-1 0 0 0;
                        1/4 -1 0 0;
                        1/4 0 -1 0;
                        0 1/4 1/4 -1]



#-----------------------------------------------------------------------------------
#Triangular - verificar se uma dada matriz é triangular
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

function is_triangular(matrix, upper = true)

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
is_triangular(upperTriangularMatrix)

#-----------------------------------------------------------------------------------
#Transposta 
#matriz tal que o elemento a_ij é igual ao elemento a_ji da matriz original
#necessária para verificar ortogonalidade

#Entrada: matriz (qualquer tamanho)
#Etapa 1: verificar se é matriz (senão retorna erro)
#Etapa 2: construir matriz transposta vazia com dimensoes corretas
#Etapa 3: popula matriz transposta onde o novo elemento a_ij é igual ao elemento a_ji da matriz original
#Saída: matriz (transposta)

function transpose(matrix)
    #checa se input é matriz
    if isa(matrix, Array{Float64,2}) == true  || isa(matrix, Array{Int64,2}) == true
        #checa se input é matriz quadrada
        numberOfLines, numberOfColumns = size(matrix)
        if numberOfLines == numberOfColumns
            #cria matriz transposta vazia com dimensoes corretas para popularmos depois
            transpose = zeros(numberOfColumns, numberOfLines)
            #o elemento a_ij da nova matriz sera o a_ji da original
            for line in 1:numberOfColumns 
                for column in 1:numberOfLines
                    transpose[line,column] = matrix[column,line]
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
    return transpose 
end

#teste
#is_triangular(transpose(upperTriangularMatrix))
transpose(upperTriangularMatrix)

#-----------------------------------------------------------------------------------
#Identidade - construir uma matriz identidade com o tamanho desejado
#identidade é propriedade de matrizes quadradas
#matriz identidade tem elementos da diagonal iguais a 1 e o resto igual a 0
#necessaria para verificar ortogonalidade

#Entrada: dimensao
#Etapa 1: construir matriz de zeros com dimensao nxn
#Etapa 2: popular elementos da diagonal com 1
#Saída: matriz  (identidade)

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

function orthogonal(matrix)
     #checa se input é matriz
    if isa(matrix, Array{Float64,2}) == true  || isa(matrix, Array{Int64,2}) == true
        #checa se input é matriz quadrada
        numberOfLines, numberOfColumns = size(matrix)
        if numberOfLines == numberOfColumns
            transposta = transpose(matrix)
            #se matrix^T * matriz = Id, é ortogonal
            if transposta*matrix == I
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

function substitution(matrix,vector) 
    #checa se input é matriz
    if isa(matrix, Array{Float64,2}) == true  || isa(matrix, Array{Int64,2}) == true
        #Se não for quadrada nxn ou vetor nao tiver dimensao correta nx1 retorna erro
        numberOfLines,numberOfColumns = size(matrix)
        if numberOfLines != numberOfColumns || numberOfLines != size(vector)[1] || 1 != size(vector)[2]
            throw(ArgumentError("'matrix' parameter should be a nxn square matrix and 
                    'vector' should be a nx1 array"))
        end
        #caso matriz seja triangular inferior
        if is_triangular(matrix,false) == true  
            #susbtituicao!
            x = zeros(numberOfLines)
            x[1] = vector[1]/matrix[1,1]
            for i in 2:numberOfLines 
            x[i] = (vector[i]- sum(matrix[i,j]*x[j] for j in 1:i-1))/matrix[i,i]
                        end 
            return x
        #Caso matriz for triangular superior
        elseif is_triangular(matrix,true) == true
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

print(A*substitution(A,b), b)

