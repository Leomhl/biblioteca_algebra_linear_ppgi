#------------------------------------------------------------------------------------------
# Á L G E B R A   L I N E A L   C O M P U T A C I O N A L
#------------------------------------------------------------------------------------------
#     J U L I A   F U N C T I O N S   L I B R A R Y
#------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
# D E S C O M P O S I Ç Ã O   D E   M A T R I Z E S
#------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
# Fatoração QR baseado em Projeções
#------------------------------------------------------------------------------------------
function Factorization_QR(A)
   m,n = size(A);
   Q = zeros(n,n);
   R = zeros(n,n);
   Q[:,1]=(1/norm(A[:,1]))*A[:,1];
   for j = 2:n
       #Obter Q por Normalizacao do vetor D
       P = zeros(n);
       for k = 1:(j-1)
           P  = P + (A[:,j]'*Q[:,k])*Q[:,k];
       end
       D = A[:,j] - P;
       Q[:,j]=(1/norm(D))*D;
   end
   R = Q'A;
   return Q, R
end
#------------------------------------------------------------------------------------------
# Fatoração L1U baseado em Aproximaçôes de Matrizes de Posto 1
# Diagonal da Matriz Triangular Inferior L con 1s
#------------------------------------------------------------------------------------------
function Factorization_L1U(A)
  m,n = size(A);
  p = min(m,n);
  L = zeros(m,p);
  U = zeros(p,n);
  T = copy(A);
  for k = 1:p
      L[:,k] = (1/T[k,k])*T[:,k];
      U[k,:] = T[k,:];
      Z = L[:,k]*U[k,:]';
      T = T - Z;
  end
  return L,U
end
#------------------------------------------------------------------------------------------
# Fatoração LU1 baseado em Aproximaçôes de Matrizes de Posto 1
# Diagonal da Matriz Triangular Superior U con 1s
#------------------------------------------------------------------------------------------
function Factorization_LU1(A)
  m,n = size(A);
  p = min(m,n);
  L = zeros(m,p);
  U = zeros(p,n);
  T = copy(A);
  for k = 1:p
      L[:,k] = T[:,k];
      U[k,:] = (1/T[k,k])*T[k,:];
      Z = L[:,k]*U[k,:]';
      T = T - Z;
  end
  return L,U
end
#----------------------------------------------------------------------
# Processo de Transformação das matrizes A e B 
# para obter uma Matriz Triangular Superior
#----------------------------------------------------------------------
function EliminacaoGaussiana_TS(A,B)
   m,n = size(A);
   for j = 1:n-1
       for i = j+1:n
           E = -1*(A[i,j]/A[j,j]);
           A[i,:] = A[i,:] + E*A[j,:];
           B[i]   = B[i]   + E*B[j];
       end
   end
   return A,B
end
#----------------------------------------------------------------------
# Processo de Transformação das matrizes A e B 
# para obter uma Matriz Triangular Inferior
#----------------------------------------------------------------------
function EliminacaoGaussiana_TI(A,B)
   m,n = size(A);
   for j = n:-1:2
       for i = 1:j-1
           E = -1*(A[i,j]/A[j,j]);
           A[i,:] = A[i,:] + E*A[j,:];
           B[i]   = B[i]   + E*B[j];
       end
   end
   return A,B
end
#----------------------------------------------------------------------
# Obtenção do vetor X de soluções do Sistema Linear AX = B
# por Substituicao Reversa usando a Matriz Triangular Superior
#----------------------------------------------------------------------
function SubstituicaoReversa_TS(A,B)
   m,n = size(A);
   X = zeros(n);
   for i = n:-1:1
       S = 0;
       for j = i+1:n
           S = S + A[i,j]*X[j];
       end
       X[i] = (B[i] - S)/A[i,i];
   end
   return X;
end
#----------------------------------------------------------------------
# Obtenção do vetor X de soluções do Sistema Linear AX = B
# por Substituicao Direta usando a Matriz Triangular Inferior
#----------------------------------------------------------------------
function SubstituicaoDireta_TI(A,B)
   m,n = size(A);
   X = zeros(n);
   for i = 1:n
       S = 0;
       for j = 1:i-1
           S = S + A[i,j]*X[j];
       end
       X[i] = (B[i] - S)/A[i,i];
   end
   return X;
end
#----------------------------------------------------------------------
# Processo de Transformação das matrizes A e B 
# para obter uma Matriz Diagonal
#----------------------------------------------------------------------
function GaussJordan(A,B)
   m,n = size(A);
   for j = 1:n
       for i = 1:n
           if i!=j
              E = -1*(A[i,j]/A[j,j]);
              A[i,:] = A[i,:] + E*A[j,:];
              B[i]   = B[i]   + E*B[j];
           end
       end
   end
   return A,B
end
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
# F U N Ç Õ E S   D E   V I S U A L I Z A Ç Ã O
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
# Formatea um Número Real em "PI" Posições Inteiras e 
# "PD" Posições Decimais
#------------------------------------------------------------------------------------------
function RealFormat(X,PI,PD)
  X = round(X; digits=PD);
  S = string(X);
  L = length(S);
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SwNegative = (S[1]=='-');
  if SwNegative==true
     if X==0
        S = "0.0";
        L = 3;
        SwNegative = false;
     else
        S = S[2:L];
        L = L - 1;
     end
  end
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  n = 0;
  Q = 0;
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  E = findfirst('e', S);
  if E==nothing
     E = 0;
  else
     Q = findfirst('-', S);
     if Q==nothing
        Q = 0;
     end
     n = parse(Int,S[(Q==0 ? E+1 : Q+1):L]);  
    #n = (Q==0 ? parse(Int,S[E+1:L]) : parse(Int,S[Q+1:L]));
     S = S[1:E-1];
     L = E-1;
  end
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  P = findfirst('.', S);
  if P==nothing
     S = string(S,".0");
     L = L + 2;
     P = L - 1;
  end
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  NI = S[1:P-1];
  LI = P - 1;
  ND = S[P+1:L];
  LD = L - P;
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if E==0
  else
     if Q==0
        ND = string(ND, "0"^(n<LD ? 0 : n-LD+1));
        LD = length(ND);
     else
        NI = string("0"^(n<LI ? 0 : n-LI+1), NI);
        LI = length(NI);
     end
     WW = NI;
     NI = (Q==0 ? string(NI,ND[1:n]) : NI[1:LI-n]);
     ND = (Q==0 ? ND[n+1:LD]         : string(WW[LI-n+1:LI],ND));
     LI = length(NI);
     LD = length(ND);
  end
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  NI = (SwNegative==true ? string("-",NI) : NI );
  LI = (SwNegative==true ? LI+1 : LI );
  NI = (PI<=LI ? NI       : string(" "^(PI-LI),NI) );
  ND = (PD<=LD ? ND[1:PD] : string(ND,"0"^(PD-LD)) );
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  return string(NI,".",ND);
end
#------------------------------------------------------------------------------------------
# Imprime Horizontalmente o vetor U considerando 
# "PI" Posições Inteiras e "PD" Posições Decimais
#------------------------------------------------------------------------------------------
function PrintVetorH(U,PI,PD)
   n = length(U)
   for k = 1:n
       print(RealFormat(U[k],PI,PD));
   end
end
#------------------------------------------------------------------------------------------
# Imprime Verticalmente o vetor U considerando 
# "PI" Posições Inteiras e "PD" Posições Decimais
#------------------------------------------------------------------------------------------
function PrintVetorV(U,PI,PD)
   n = length(U)
   for k = 1:n
       println(RealFormat(U[k],PI,PD));
   end
end
#------------------------------------------------------------------------------------------
# Imprime os elementos da matriz A considerando 
# "PI" Posições Inteiras e "PD" Posições Decimais
#------------------------------------------------------------------------------------------
function PrintMatrix(A,PI,PD)
   m,n = size(A)
   for i = 1:m
       for j = 1:n
           print(RealFormat(A[i,j],PI,PD));
       end
       println("");
   end
end
#------------------------------------------------------------------------------------------
# Imprime os elementos da matriz A e do vetor B 
# considerando "PI" Posições Inteiras e "PD" Posições Decimais
#------------------------------------------------------------------------------------------
function PrintLinearSystem(A,B,PI,PD)
   m,n = size(A)
   for i = 1:m
       for j = 1:n
           print(RealFormat(A[i,j],PI,PD));
       end
       print(" "^PI,"│");
       print(RealFormat(B[i],PI,PD));
       println("");
   end
end
#------------------------------------------------------------------------------------------
# Imprime a Matriz em um Arquivo TXT
#------------------------------------------------------------------------------------------
function MatrixToText(A,PI,PD,FILE)
   m,n = size(A);
   F = open(FILE,"w");
   for i = 1:m
       for j = 1:n
           write(F,RealFormat(A[i,j],PI,PD));
       end
       write(F,"\n");
   end
   close(F);
end


