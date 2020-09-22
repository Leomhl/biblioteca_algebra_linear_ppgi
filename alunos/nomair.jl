#    ÁLGEBRA LINEAR COMPUTACIONAL
#
#    BIBLIOTECA DE FUNÇÕES DO JULIA
#
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

# Funções Auxiliares

#1 Se v for um matriz de 1 coluna/linha, transforma v em um vetor
function ensure_vector(v) 
    @assert ndims(v) == 1 || (max_val(size(v)) == length(v))
    return v[:]
end

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

#4 Gaussian Elimunations usa uma sequência de operações de linha elementares para modificar a matriz até que o canto esquerdo 
# inferior da matriz seja preenchido com zeros, tanto quanto possível

function gaussian_eliminations(A)
    rows = size(A,1)
    cols = size(A,2)
    fraction = A[1,1]/A[1,1]
    # Row index
    row = 1
    U = copy!(similar(A, typeof(fraction)), A)
    # Main loop going through all columns
    for col = 1:(cols-1)
        
        # finding the maximum element for each column
        max_index = argmax(abs.(A[row:end,col])) + row-1
        # Check to make sure matrix is good!
        if (U[max_index, col] == 0)
            continue
        end

        # swap row with highest value for that column to the top
        temp_vector = U[max_index, :]
        U[max_index, :] = U[row, :]
        U[row, :] = temp_vector
        
        # Loop for all remaining rows
        for i = (row+1):rows
            # finding fraction
            fraction = U[i,col]/U[row,col]
            # loop through all columns for that row
            for j = (col+1):cols
                 # re-evaluate each element  
                 U[i,j] -= U[row,j]*fraction
            end
            # Set lower elements to 0
            U[i,col] = 0
        end
        row += 1
    end
    return U
end

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------


# Funções Principais

#------------------------------------------------------------------------------------------
#1 Se 2 vetores são perpendiculares
# u*v = 0
#------------------------------------------------------------------------------------------

function vetores_perpendiculares(u,v)
    result = dot_prod(u,v)
    if (result == 0) 
        return true
    else
        return false
    end
end

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------
#2 Ângulo entre dois vetores
# cos theta = <u,v>/norm(u)norm(v)
#------------------------------------------------------------------------------------------
function vetores_angulo(u,v)
    u1 = ensure_vector(u)
    v1 = ensure_vector(v)
    result = 0
    result = dot_prod(u,v)/(norm(u1)*norm(v1))
    return result
end

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------
#3 Posto da Matriz
# Corresponde ao número de linhas ou colunas linearmente independentes da matriz.
#------------------------------------------------------------------------------------------

function matrix_rank(A)
    U = gaussian_eliminations(A)
    rows = size(U,1)
    rank = rows
    for i = 1:rows
        row_sum =0
        for j = 1:length(U[i, :])
            row_sum += U[i,j]
        end
        if (row_sum ==0) 
            rank-=1
        end
    end
    return rank
end

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

#testes pra ver se as funções estão corretas
function testes()
    # using LinearAlgebra
    tol = 10e-5;
    
    A = rand(10, 10);
    @assert matrix_rank(A) - rank(A) < tol
    

    u = rand(10);
    v = rand(10);
    @assert dot(u, v)/norm(u)norm(v) - vetores_angulo(u,v) < tol
    
    
    r = [ 2 , 1, 4,2];
    s = [ 3 , -6, -1,2];
    @assert vetores_perpendiculares(r,s) == true

    r = [ 2 , 1, 4,2];
    s = [ 3 , -6, -1,-2];
    @assert vetores_perpendiculares(r,s) == false
    
end

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------


