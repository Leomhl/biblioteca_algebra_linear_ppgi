#******************************************************************************************
#Funções a serem implementadas:
#Triangular - verificar se uma dada matriz é triangular
#Ortogonal - verificar se uma dada matriz é ortogonal
#Ortonormal - verificar se um dado conjunto de vetores forma base ortonormal
#Substituição reversa - para solução de sistemas de equações lineares
#Pseudo inversa - utilizando decomposição svd
#SVD - decomposição (para poder fazer pseudo inversa)
#------------------------------------------------------------------------------------------

#******************************************************************************************
#MATRIZ TRIANGULAR - verificar se uma dada matriz é triangular
#é chamada triangular inferior se todas as entradas acima da diagonal principal são 0.
#é chamada triangular superior se todas as entradas abaixo da diagonal principal são 0

#Entrada: matriz, se teste é para superior/inferior (default é superior, para inferior usar false em upper)
#Etapa 1: transforma entrada em array (caso seja de outro tipo); 
#Etapa 2: olha elementos abaixo ou acima da diagonal principal para determinar se é triangular ou não
#Saída: True (é triangular) ou False ( não é triangular)
#------------------------------------------------------------------------------------------

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



#******************************************************************************************
# MATRIZ TRANSPOSTA
#matriz tal que o elemento a_ij é igual ao elemento a_ji da matriz original
#necessária para verificar ortogonalidade

#Entrada: matriz (qualquer tamanho)
#Etapa 1: transforma entrada em array
#Etapa 2: construir matriz transposta vazia com dimensoes corretas
#Etapa 3: popula matriz transposta onde o novo elemento a_ij é igual ao elemento a_ji da matriz original
#Saída: matriz (transposta)
#------------------------------------------------------------------------------------------

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

#******************************************************************************************
# MATRIZ IDENTIDADE
# Retorna a identidade de uma matriz (digianonal pricipal = 1, restante = 0)

# Entrada: tamanho da matriz a ser construida
# Saídas: matriz identidade
# Autor: Gabriella
#------------------------------------------------------------------------------------------
function identidade(tamanho)
    #constroi matriz de 0 do tamanho correto
    I=zeros(tamanho,tamanho) 
    #preenche diagonal com 1
    for i=1:tamanho
        I[i,i]=1 
    end 
    return I      
end

#******************************************************************************************
#MATRIZ ORTOGONAL - verificar se uma matriz é ortogonal
#Note que ser ortogonal é uma propriedade de matrizes quadradas
#portanto, a matriz do input deverá ser quadrada
#Uma matriz é ortogonal quando suas linhas e colunas são vetores ortogonais
#outra forma de saber é se a transposta vezes a matriz é igual identidade

#Entrada: matriz
#Etapa 1: verifica se entrada é matriz quadrada, senao retorna erro
#Etapa 2: calcula transposta da matriz
#Etapa 3: multiplica transposta pela matriz original e verifica se deu identidade
#Saída: true(é ortogonal) ou false(não é ortogonal)
#------------------------------------------------------------------------------------------

function ortogonal(matriz)
    #transforma entrada em array
    matriz = Array(matriz)
    #checa se input é matriz quadrada
    numeroDeLinhas, numeroDeColunas = size(matriz)
    if numeroDeLinhas == numeroDeColunas
        matrizTransposta = transposta(matriz)
        #se matrix^T * matriz = Id, é ortogonal
        if matrizTransposta*matriz == identidade(numeroDeLinhas)
            return true
        else
            return false
        end
    else
        #retorna erro se matriz nao e quadrada
        throw(ArgumentError("Matriz deve ser quadrada. Matrizes retangulares podem ter colunas OU linhas ortonormais.
                Porém não os dois ao mesmo tempo."))
    end
end


#******************************************************************************************
#BASE ORTONORMAL - verificar se um dado conjunto de vetores forma base ortonormal
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
#------------------------------------------------------------------------------------------

function base_ortonormal(conjunto_de_vetores)
    #transforma entrada em array
    conjunto_de_vetores = Array(conjunto_de_vetores)
        #pega tamanho do primeiro vetor, iremos comparar com o resto
        numeroDeLinhas, numeroDeColunas = size(conjunto_de_vetores[1])
        #alguma das dimensoes precisa ser 1, senao nao é vetor...
        if numeroDeLinhas != 1 && numeroDeColunas != 1
            #se nenhuma for 1, retorna erro
            throw(ArgumentError("Elementos do conjunto de vetores devem ser vetores (dimensao nx1 ou 1xn)"))
        else
            #para cada vetor, ele deve ter o mesmo tamanho do primeiro e comprimento 1
            for vetor in conjunto_de_vetores
                if numeroDeLinhas == size(vetor)[1] && numeroDeColunas == size(vetor)[2] && 
                    norma(vetor) <= 1.001 && norma(vetor) >= 0.999
                    continue
                else
                    return false
                #cada par de vetor devera ser ortogonal entre si
                for vetor2 in conjunto_de_vetores
                    if produto_interno(vetor,vetor2) == 0   
                        continue
                    else
                        return false     
                    end
                end
                end
            end
        end
    #se chegamos ate aqui, é pq a base é ortonormal
    return true
end

#******************************************************************************************
#SUBSTITUICAO - solucionar sistemas de equaçoes lineares Ax=b 
#a matriz devera ser quadrada e triangular

#Entrada: matriz  A e vetor b
#Etapa 1: verificar se matriz QUADRAD é triangular superior ou inferior
#Etapa 2: verificar se vetor tem dimensao correta 
#Etapa 3: substituiçao
#Saida: vetor x - solucao do sistema linear
#------------------------------------------------------------------------------------------

function substituicao(matriz,vetor) 
    #transforma matriz em array
    matriz = Array(matriz)
    #Se não for quadrada nxn ou vetor nao tiver dimensao correta nx1 retorna erro
    numeroDeLinhas,numeroDeColunas = size(matriz)
    if numeroDeLinhas != numeroDeColunas || numeroDeLinhas != size(vetor)[1] || 1 != size(vetor)[2]
        throw(ArgumentError("Matriz deve ser quadrada nxn e vetor deve ser nx1"))
    end
    #caso matriz seja triangular inferior
    if e_triangular(matriz,false) == true  
        #susbtituicao!
        x = zeros(numeroDeColunas)
        x[1] = vetor[1]/matriz[1,1]
        for i in 2:numeroDeLinhas 
        x[i] = (vetor[i]- sum(matriz[i,j]*x[j] for j in 1:i-1))/matriz[i,i]
                    end 
        return x
    #Caso matriz for triangular superior
    elseif e_triangular(matriz,true) == true
        #substituicao!
        x = zeros(numeroDeColunas) 
        x[numeroDeLinhas] = vetor[numeroDeLinhas]/matriz[numeroDeLinhas,numeroDeColunas]
        for i in numeroDeLinhas-1:-1:1
            x[i] = (vetor[i]- sum(matriz[i,j]*x[j] for j in i+1:numeroDeLinhas))/matriz[i,i]
                    end            
        return x
    else
        throw(ArgumentError("Matriz deve ser triangular"))
    end
end 
 

#******************************************************************************************
#TESTE TRIANGULAR
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

matriz_gastao = lowerTriangularMatrix'

e_triangular(matriz_gastao)

e_triangular(squareTriL,false)

e_triangular(squareTriL)

e_triangular(rectangularTriU)

e_triangular(rectangularTriL, false)

e_triangular(rectangularTriL)

e_triangular(notTriangularMatrix, Array(matriz))

#******************************************************************************************
#TESTE TRANSPOSTA - se for triangular inferior deve se tornar triangular superior!
e_triangular(transposta(squareTriL))


transposta(squareTriU)

#******************************************************************************************
#TESTE ORTOGONAL
ortogonal([1 0; 0 1])

ortogonal([1 1; 0 0])

#******************************************************************************************
#TESTE BASE ORTONORMAL
conjunto_de_vetores_ortonormal = [[1/sqrt(2) 0 -1/sqrt(2)], [1/2 sqrt(2)/2 1/2], [1/2 -sqrt(2)/2 1/2]]
conjunto_de_vetores_ortogonal = [[1 1 1], [-2 1 1], [0 1 -1]]

base_ortonormal(conjunto_de_vetores_ortogonal)

base_ortonormal(conjunto_de_vetores_ortonormal)

#******************************************************************************************
#TESTE SUBSTITUICAO
b = randn(3,1)
A = [1.0      0.0      0.0;
 2.12374  1.0      0.0;
 1.20021  1.04731  1.0]

#tem que dar a mesma coisa
print(A*substituicao(A,b), b)


b = randn(3,1)
A = [1.0      2.0      0.0;
 2.12374  1.0      0.0;
 1.20021  1.04731  1.0]

#é para dar erro pois nao é triangular
print(A*substituicao(A,b), b)

b = randn(4,1)

print(squareTriU*substituicao(squareTriU,b),b)


b = randn(4,1)

print(rectangularTriU*substituicao(rectangularTriU,b),b)

b = randn(4,1)

print(squareTriL*substituicao(squareTriL,b),b)
