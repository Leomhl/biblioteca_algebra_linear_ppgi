#    ÁLGEBRA LINEAR COMPUTACIONAL
#
#    TESTES PARA A BIBLIOTECA DE FUNÇÕES DO JULIA
#
#-------------------------------------------------------------------------------
#
#    DECOMPOSIÇÃO QR
#
function teste_qr(k)
    Tudo_certo=true
    for n=1:k
        for m=1:k
            A=randn(n,m)
            for i=1:2
                Q,R=decomposição_qr(A)
                if norm(A-Q*R)>0.00001
                    Tudo_certo=false
                end
                A=rand(m,n)
            end
        end
    end
        
    return Tudo_certo
end
#
#-------------------------------------------------------------------------------
#
#    DECOMPOSIÇÃO LU
#
function teste_lu(k)
    Tudo_certo=true
    for n=1:k
        for m=1:k
            A=randn(n,m)
            for i=1:2
                L,U=decomposição_lu(A)
                if norm(A-L*U)>0.00001
                    Tudo_certo=false
                end
                A=rand(m,n)
            end
        end
    end
    
    return Tudo_certo
end
#
#-------------------------------------------------------------------------------
#
#    Resolução de sistema linear
#
function test()
    n=Int(abs(rand(Int8)))
    A=rand(n,n)
    b=rand(n)    
    X=resolver(A,b)
    ok=norm(A*X-b)
    return X,ok
end