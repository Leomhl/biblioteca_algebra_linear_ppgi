#    ÁLGEBRA LINEAR COMPUTACIONAL
#
#    BIBLIOTECA DE FUNÇÕES DO JULIA
#
#
#-------------------------------------------------------------------------------
#
#    FUNÇÕES BÁSICAS
#
function somatorio(v)
    u=ones(length(v))
    return v'u
end
#
function media(v)
    return (1/length(v))*somatório(v)
end
#
function media_ponderada(v,w)
    if length(v)!=length(w)
        print("A quantidade de dados deve ser a mesma quantidade de pesos.")
    else
        return (1/somatório(w))*v'w
    end
end
#
function produto_interno(u,v)
    if length(v)!=length(w)
        print("As quantidades de coordenadas dos dois vetores devem ser as mesmas!!!")
    else
        return u'v
    end
end
#
angulo(x,y)=acos(x'*y/(norm(x)*norm(y)))
#
function norma(v)
    return sqrt(v'v)
end
#
function normalizar(v)
    return v/norma(v)
end
#
function distancia(u,v)
    return norma(u-v)
end
#
function vetor_canonico(i,n)
    if 0<i<=n
        return [zeros(i-1);1;zeros(n-i)]
        else println("Valores inválidos.") 
    end
end
#
function combinacao_linear(escalares,vetores)
    return escalares.*vetores
end
#
#-------------------------------------------------------------------------------
#
#    PROJEÇÃO
#
function projeção(a,v)
    v=v/norm(v)
    return (a'v)*v
end
#
#-------------------------------------------------------------------------------
#
#    DECOMPOSIÇÃO QR
#
function decomposição_qr(A)
    n,m=size(A)
    dim=min(n,m)
    Q=zeros(n,dim)
    R=zeros(dim,m)
    for i=1:dim
        
        Q[:,i]=A[:,i]/norm(A[:,i])
        
        R[i,:]=Q[:,i]'*A
        
        S=Q[:,i]*R[i,:]'    #S é a matriz a ser subtraída de A

        A=A-S        
    end
    
   return Q,R 
end
#
function decomposição_qr2(A)
    n,m=size(A)
    dim=min(n,m)
    Q=zeros(n,dim)
    R=zeros()
    for i=1:dim
        
        Q[:,i]=A[:,1]/norm(A[:,i])
        R[i,:]=Q[:,i]'*A
        
    end
    
   return Q,R 
end
#
#-------------------------------------------------------------------------------
#
#    DECOMPOSIÇÃO LU
#
function decomposição_lu(A)
    n,m=size(A)
    dim=min(n,m)
    L=zeros(n,dim)
    U=zeros(dim,m)
    
    for i=1:dim
        L[:,i]=A[:,i]/A[i,i]       #Coluna i de L

        U[i,:]=A[i,:]       #Linha i de U
        
        S=L[:,i]*A[i,:]'    #S é a matriz que iremos subtrair: L[:,i]*U[i,:]

        A=A-S
    end
    
    return L,U
end
#
function decomposição_lu(A)
    Ad=copy(A)              #A decomposta
    n=length(A[1,:])
    L=zeros(n,n)
    U=zeros(n,n)
    
    for i=1:n
        L[:,i]=Ad[:,i]       #linha i de L

        D=L[:,i]*Ad[i,:]'    #D é a matrir L[:,i]*U[i,:]

        U[i,:]=D[i,:]        #Coluna i de U

        Ad=Ad-D
    end
    U[n,:]=[zeros(n-1);1]
    return L,U
end
#
#-------------------------------------------------------------------------------
#
#    SUBSTITUIÇÃO REVERSA
#
function substituição_triangular_inferior(A,b)
    n=length(b)
    x=zeros(n)
    D=0
    x[1]=b[1]/A[1,1]
    for i=1:n
        D=A[i,1:(i-1)]'*x[1:(i-1)]
        x[i]=(b[i]-D)/A[i,i]
    end
        
    return x
end
#
function subs_reversa(A,b)
 
    n=length(b)
    x=zeros(n,1)
   
    for k=n:1 #linhas
        x[k]=b[k]/A[k,k]
   
        x[k]=(b[k]-A[2,3]*x[3])/A[k,k]
   
        x[k]=(b[k]-A[1,2]*x[2]-A[1,3]*x[3])/A[k,k]
    end
    return x
end   
#
function subs_reversa2(A,b)
  
    n=length(b)
    x=zeros(n,1)
    
    for k=n:1 #linhas
        x[k]=b[k]/A[k,k]
    
        x[k]=(b[k]-A[2,3]*x[3])/A[k,k]
    
        x[k]=(b[k]-A[1,2]*x[2]-A[1,3]*x[3])/A[k,k]
    end
    return x
end    
#-------------------------------------------------------------------------------
#
#    SUBSTITUIÇÃO NORMAL
#
function substituição_triangular_superior(A,b)
    n=length(b)
    x=zeros(n)
    D=0
    x[n]=b[n]/A[n,n]
    for i=n:-1:1
        D=A[i,1:(i-1)]'*x[1:(i-1)]
        x[i]=(b[i]-D)/A[i,i]
    end
    
    return x 
end
#
#-------------------------------------------------------------------------------
#
#    RESOLVER SISTEMA LINEAR POR DECOMPOSIÇÃO LU
#

function resolver_LU(A,b)
    
    L,U=decomposição_lu(A)
    y=substituição_triangular_inferior(L,b)
    x=substituição_triangular_superior(U,y)
    
    return x
end
#
#-------------------------------------------------------------------------------
#
#    ELIMINAÇÃO GAUSSIANA
#
function eliminacao(A,b)
    n=length(b)
   
    for i=1:(n-1)
   
        for k = (i+1):n # na coluna
            b[k]=b[k]-(A[k,i]/A[i,i])*b[i]
            A[k,:]=A[k,:]-(A[k,i]/A[i,i])*A[i,:]
        end  
   
    end    
    return A,b
end  
#
function eliminacao2(A,b)
    n=length(b)
    
    for i=1:(n-1)
    
        for k = (i+1):n # na coluna
            b[k]=b[k]-(A[k,i]/A[i,i])*b[i]
            A[k,:]=A[k,:]-(A[k,i]/A[i,i])*A[i,:]
        end  
    
    end    
    return A,b
end   
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#    SVD
#-------------------------------------------------------------------------------
function erro(a,v)
    p=projeção(a,v)
    return norm(a-p)
end

function erro_total(pontos,v)
    erro=0
    for i=1:length(pontos[:,1])
        p=projeção(pontos[i,:],v)
        erro+=erro(p-v)
    end
    return erro
end

function matriz_erro(pontos,v)
    return pontos'-projeção(pontos,v)
end
#-------------------------------------------------------------------------------
# Preciso das funções: norma; multiplicação matricial; transposta.
#Input: uma matriz com os pontos dados
#Output: o vetor direção "v" da reta que melhor se aproxima dos pontos dados
function SVD(pontos)
    A=copy(pontos) #matriz criada com os pontos dados (Tem que testar se precisa pegar a transposta.)
    v=A[:,1]       
    v=v/norm(v)    #Primeiro candidato a "v" é o primeiro ponto (normalizado) dado (vendo os pontos como os vetores colunas de A).
    erro=0         
    erro_novo=norm(A-v*(A'*v)')  #O erro é verificado vendo a "distância" entre os pontos dados e a reta dada pelo vetor v".
    
    println(erro_novo)
    
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
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#
#    SISTEMAS DINÂMICOS
#
function dinamica(A,x0,k)
    x=x_0
    for i=1:k
        x=A*x
    end

    return x
end
