#    ÁLGEBRA LINEAR COMPUTACIONAL
#
#    TESTES PARA A BIBLIOTECA DE FUNÇÕES DO JULIA




#--------------------------------------------------------------------------------------------------------------------------
#    TESTES COM MATRIZES ALEATÓRIAS
#------------------------------------------------------------------------------------------
# A função faz um teste com matrizes variando o número de linhas e colunas até o valor máximo de k.
# Entradas: número inteiro positivo, função 
# Saída: booleano
#------------------------------------------------------------------------------------------
#function testes_matrizes(n,m,teste)	# "teste" é a função teste que será usada.
#    Tudo_certo=true
#    for n=1:k
#        for m=1:k
#            A=randn(n,m)
#            if teste(A)
#                Tudo_certo=false
#            end
#        end
#    end
#    
#    return Tudo_certo
#end
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------




#--------------------------------------------------------------------------------------------------------------------------
#    DECOMPOSIÇÃO QR
#------------------------------------------------------------------------------------------

#
#
#------------------------------------------------------------------------------------------
function teste_qr(k)
    Tudo_certo=true
    for n=1:k
        for m=1:k
            A=randn(n,m)
            Q,R=decomposição_qr(A)
            if norm(A-Q*R)>0.00001
                if Q*Q'!=I
                    Tudo_certo=false
                end
            end
        end
    end
        
    return Tudo_certo
end
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------




#--------------------------------------------------------------------------------------------------------------------------
#    DECOMPOSIÇÃO LU
#------------------------------------------------------------------------------------------
#
#
#
#------------------------------------------------------------------------------------------
function teste_lu(k)
    Tudo_certo=true
    for n=1:k
        for m=1:k
            A=randn(n,m)
            L,U=decomposição_lu(A)
            if norm(A-L*U)>0.00001
                # L é triangular inferior?
                # U é triangular superior?
                Tudo_certo=false
            end
        end
    end
    
    return Tudo_certo
end
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------




#--------------------------------------------------------------------------------------------------------------------------
#    Resolução de sistema linear
#------------------------------------------------------------------------------------------
#
#    
#
function test()
    n=Int(abs(rand(Int8)))
    A=rand(n,n)
    b=rand(n)    
    X=resolver(A,b)
    ok=norm(A*X-b)
    return X,ok
end
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------




#--------------------------------------------------------------------------------------------------------------------------
#    SVD
#------------------------------------------------------------------------------------------
function teste_SVD()
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
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------




#--------------------------------------------------------------------------------------------------------------------------
#    SISTEMAS DINÂMICOS
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
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------





