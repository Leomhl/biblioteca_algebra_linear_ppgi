#------------------------------------------------------------------------------------------
# Á L G E B R A   L I N E A L   C O M P U T A C I O N A L
#------------------------------------------------------------------------------------------
#     J U L I A   F U N C T I O N S   L I B R A R Y
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
function NormVector(V)
   n = length(V);
   S = 0;
   for i = 1:n
       S = S + (V[i]^2);
   end
   return sqrt(S);
end
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
function NormMatrix(A)
   m,n = size(A);
   S = 0;
   for i = 1:m
       for j = 1:n
           S = S + (A[i,j]^2);
       end
   end
   return sqrt(S);
end
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
function LU(A)
  return Factorization_L1U(A);
end
function LU(A,Sw)
  if Sw
     return Factorization_L1U(A);
  else
     return Factorization_LU1(A);
  end 
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
   X = zeros(n);
   for j = 1:n
       for i = 1:n
           if i!=j
              E = -1*(A[i,j]/A[j,j]);
              A[i,:] = A[i,:] + E*A[j,:];
              B[i]   = B[i]   + E*B[j];
           end
       end
       X[j] = B[j]/A[j,j];
   end
   return X;
end
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
# Metodo de Jacobi
# Realiza Aproximações Sucesivas para obter a Solução
# Usa o vetor auxiliar Y para guardar a solucao da iteracao anterior
#------------------------------------------------------------------------------------------
function Jacobi(A,B,E)
   m,n = size(A);
   Y = rand(n);
   X = zeros(n);
   t = 0;
   while NormVector(X-Y)>E
      t = t + 1;
      Y = copy(X)
      for i = 1:n
          S = 0;
          for j = 1:n
              if i!=j
                 S = S + A[i,j]*Y[j];
              end
          end
          X[i] = (B[i]- S)/A[i,i];
      end
   end
   return X
end
#------------------------------------------------------------------------------------------
# Metodo de Gauss-Seidel
# Realiza Aproximações Sucesivas para obter a Solução
# Usa o vetor auxiliar Y para guardar a solucao da iteracao anterior
#------------------------------------------------------------------------------------------
function Gauss_Seidel(A,B,E)
   m,n = size(A);
   Y = rand(n);
   X = zeros(n);
   t = 0;
   while NormVector(X-Y)>E
      t = t + 1;
      Y = copy(X)
      for i = 1:n
          S = 0;
          for j = 1:n
              if i!=j
                 S = S + A[i,j]*X[j];
              end
          end
          X[i] = (B[i]- S)/A[i,i];
      end
   end
   return X
end
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
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Gera Codigo CSS para formato da Pagina HTML
#-------------------------------------------------------------------------------
function GetCSS()
   ForeColorCell     = "#000000";  #Text Color of Cell
   BackColorCell     = "#FFDB4D";  #Background Color of Cell
   ForeColorIndex    = "#FFFFFF";  #Text Color of Index
   BackColorIndex    = "#505050";  #Background Color of Index
   ForeColorDiagonal = "#000000";  #Text Color of Diagonal
   BackColorDiagonal = "#99E6FF";  #Background Color of Diagonal
   TitleColor        = "#0000FF";  #Title Color
   TextColor         = "#00FF00";  #Text Color
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   S = "<STYLE>\n";
   S = string(S,".TABLE     {  padding: 10px; border-spacing:1px 1px; }\n" );
   S = string(S,".INDEX     {  color: " , ForeColorIndex    , "; background-color: " , BackColorIndex    , "; text-align: center;  padding: 10px;  font-family: Calibri, Tahoma, Verdana; }\n" );
   S = string(S,".DIAGONAL  {  color: " , ForeColorDiagonal , "; background-color: " , BackColorDiagonal , "; text-align: center;  padding: 15px;  font-family: 'Courier New', Tahoma, Verdana; font-weight: bold; }\n" );
   S = string(S,".CELL      {  color: " , ForeColorCell     , "; background-color: " , BackColorCell     , "; text-align: center;  padding: 15px;  font-family: 'Courier New', Tahoma, Verdana; }\n" );
   S = string(S,".TITLE     {  color: " , TitleColor        , "; text-align: left;  padding: 15px;  font-family: Tahoma, Verdana, Calibri;  font-size: 28px; }\n" );
   S = string(S,".TEXT      {  color: " , TextColor         , "; text-align: left;  padding: 15px;  font-family: Tahoma, Verdana, Calibri;  font-size: 28px; }\n</STYLE>\n" );
   S = string(S,"</STYLE>\n" );
   #-------------------------------------------------------------------------------
   return string(S);
end
#-------------------------------------------------------------------------------
# Gera o Codigo HTML associado para a Matriz 
#-------------------------------------------------------------------------------
function WriteMatrix(F,A,PD,ShowIndex,ShowDiagonal)
   m,n = size(A);
   #-------------------------------------------------------------------------------
   write(F,"<TABLE class=TABLE>\n");
   if ShowIndex==true
      write(F,"<TR><TD> </TD>");
      for j = 1:n
          write(F,"<TD CLASS=INDEX>",string(j),"</TD>" );
      end
      write(F,"</TR>\n");
   end
   for i = 1:m
       write(F,"<TR ALIGN=CENTER>");
       if ShowIndex==true
          write(F,"<TD CLASS=INDEX>",string(i),"</TD>" );
       end
       for j = 1:n
           write(F, (ShowDiagonal==true && i==j ? "<TD CLASS=DIAGONAL>" : "<TD CLASS=CELL>"), RealFormat(A[i,j],0,PD), "</TD>" );
       end
       write(F,"</TR>\n");
   end
   write(F,"</TABLE><br><br>\n</BODY>\n</HTML>\n");
end
#-------------------------------------------------------------------------------
# Resultados da Matriz em Arquivo HTML
#-------------------------------------------------------------------------------
function MatrixToHTML(A,PD,FILE,ShowIndex,ShowDiagonal,Title)
   #-------------------------------------------------------------------------------
   F = open(FILE,"w");
   #-------------------------------------------------------------------------------
   write(F,"<HTML>\n<HEAD>\n");
   write(F,GetCSS());
   write(F,"</HEAD>\n<BODY>\n");
   #-------------------------------------------------------------------------------
   write(F,"<H1>",Title,"</H1>\n","<HR>\n");
   #-------------------------------------------------------------------------------
   WriteMatrix(F,A,PD,ShowIndex,ShowDiagonal);
   #-------------------------------------------------------------------------------
   close(F);
end
#-------------------------------------------------------------------------------
# Simulação da Fatoração QR
#-------------------------------------------------------------------------------
function Simulation_QR(A,PD,FILE)
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   F = open(FILE,"w");
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write(F,"<HTML>\n<HEAD>\n");
   write(F,GetCSS());
   write(F,"</HEAD>\n<BODY>\n");
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write(F,"<H1>Simula&ccedil;&atilde;o do M&eacute;todo QR</H1>\n" );
   write(F,"<HR>\n");
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write(F,"<SPAN CLASS=TITLE>Matriz A</SPAN>\n");
   WriteMatrix(F,A,PD,true,false);
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   m,n = size(A);
   Q = zeros(n,n);
   R = zeros(n,n);
   Q[:,1]=(1/norm(A[:,1]))*A[:,1];
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write(F,"<SPAN CLASS=TITLE>Passo 1: Inicio</SPAN>\n<HR>\n");
   write(F,"<SPAN CLASS=TITLE>Matriz Q</SPAN>\n");
   WriteMatrix(F,Q,PD,true,false);
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   for j = 2:n
       P = zeros(n)
       for k = 1:(j-1)
           P  = P + (A[:,j]'*Q[:,k])*Q[:,k]
       end         
       D = A[:,j] - P
       Q[:,j]=(1/norm(D))*D
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       write(F,"<SPAN CLASS=TITLE>Passo ",string(j),"</SPAN>\n<HR>\n");
       write(F,"<SPAN CLASS=TITLE>Matriz Q</SPAN>\n");
       WriteMatrix(F,Q,PD,true,false);
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   end
   R = Q'A 
   write(F,"</BODY>\n</HTML>\n");
   close(F);
end
#-------------------------------------------------------------------------------
# Simulação da Fatoração L1U por Aproximação de Matriz de Posto 1
#-------------------------------------------------------------------------------
function Simulation_L1U(A,PD,FILE)
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   F = open(FILE,"w");
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write(F,"<HTML>\n<HEAD>\n");
   write(F,GetCSS());
   write(F,"</HEAD>\n<BODY>\n");
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write(F,"<H1>Simula&ccedil;&atilde;o do M&eacute;todo L1U por Aproxima&ccedil;&atilde;o de Matriz de Posto 1</H1>\n" );
   write(F,"<HR>\n");
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write(F,"<SPAN CLASS=TITLE>Matriz A</SPAN>\n");
   WriteMatrix(F,A,PD,true,true);
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   m,n = size(A);
   p = min(m,n);
   L = zeros(m,p);
   U = zeros(p,n);
   T = copy(A);
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write(F,"<SPAN CLASS=TITLE>Passo 0: Inicio</SPAN>\n<HR>\n");
   write(F,"<SPAN CLASS=TITLE>Matriz L</SPAN>\n");
   WriteMatrix(F,L,PD,true,true);
   write(F,"<SPAN CLASS=TITLE>Matriz U</SPAN>\n");
   WriteMatrix(F,U,PD,true,true);
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   for k = 1:p
       L[:,k] = (1/T[k,k])*T[:,k];
       U[k,:] = T[k,:];
       Z = L[:,k]*U[k,:]';
       T = T - Z;
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       write(F,"<SPAN CLASS=TITLE>Passo ",string(k),"</SPAN>\n<HR>\n");
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       write(F,"<SPAN CLASS=TITLE>Matriz de Posto 1</SPAN>\n");
       WriteMatrix(F,Z,PD,true,true);
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       write(F,"<SPAN CLASS=TITLE>Matriz Resultante</SPAN>\n");
       WriteMatrix(F,T,PD,true,true);
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       write(F,"<SPAN CLASS=TITLE>Matriz L</SPAN>\n");
       WriteMatrix(F,L,PD,true,true);
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       write(F,"<SPAN CLASS=TITLE>Matriz U</SPAN>\n");
       WriteMatrix(F,U,PD,true,true);
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   end
   write(F,"</BODY>\n</HTML>\n");
   close(F);
end
#-------------------------------------------------------------------------------
# Simulação da Fatoração LU1 por Aproximação de Matriz de Posto 1
#-------------------------------------------------------------------------------
function Simulation_LU1(A,PD,FILE)
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   F = open(FILE,"w");
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write(F,"<HTML>\n<HEAD>\n");
   write(F,GetCSS());
   write(F,"</HEAD>\n<BODY>\n");
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write(F,"<H1>Simula&ccedil;&atilde;o do M&eacute;todo LU1 por Aproxima&ccedil;&atilde;o de Matriz de Posto 1</H1>\n" );
   write(F,"<HR>\n");
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write(F,"<SPAN CLASS=TITLE>Matriz A</SPAN>\n");
   WriteMatrix(F,A,PD,true,true);
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   m,n = size(A);
   p = min(m,n);
   L = zeros(m,p);
   U = zeros(p,n);
   T = copy(A);
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write(F,"<SPAN CLASS=TITLE>Passo 0: Inicio</SPAN>\n<HR>\n");
   write(F,"<SPAN CLASS=TITLE>Matriz L</SPAN>\n");
   WriteMatrix(F,L,PD,true,true);
   write(F,"<SPAN CLASS=TITLE>Matriz U</SPAN>\n");
   WriteMatrix(F,U,PD,true,true);
   #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   for k = 1:p
       L[:,k] = T[:,k];
       U[k,:] = (1/T[k,k])*T[k,:];
       Z = L[:,k]*U[k,:]';
       T = T - Z;
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       write(F,"<SPAN CLASS=TITLE>Passo ",string(k),"</SPAN>\n<HR>\n");
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       write(F,"<SPAN CLASS=TITLE>Matriz de Posto 1</SPAN>\n");
       WriteMatrix(F,Z,PD,true,true);
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       write(F,"<SPAN CLASS=TITLE>Matriz Resultante</SPAN>\n");
       WriteMatrix(F,T,PD,true,true);
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       write(F,"<SPAN CLASS=TITLE>Matriz L</SPAN>\n");
       WriteMatrix(F,L,PD,true,true);
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       write(F,"<SPAN CLASS=TITLE>Matriz U</SPAN>\n");
       WriteMatrix(F,U,PD,true,true);
       #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   end
   write(F,"</BODY>\n</HTML>\n");
   close(F);
end
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#----------------------------------------------------------------------
# P R O C E S S O S   D E   T E S T E   E   S I M U L A Ç Ã O
#----------------------------------------------------------------------
# Simulação da EliminacaoGaussiana baseado em Matriz Triangular Superior
#----------------------------------------------------------------------
function Teste_EliminacaoGaussiana_TS(AA,BB,PI,PD)
   println("---------------------------------------------------------------------------");
   println("T R A N S F O R M A Ç Ã O   P O R   E L I M I N A Ç Ã O   G A U S S I A N A");
   println("---------------------------------------------------------------------------");
   println("\nDado o Sistema Linear Ax=B\n");
   A = copy(AA);
   B = copy(BB);
   PrintLinearSystem(A,B,PI,PD);
   m,n = size(A);
   X = zeros(n);
   for j = 1:n-1
       for i = j+1:n
           E = -1*(A[i,j]/A[j,j]);
           A[i,:] = A[i,:] + E*A[j,:];
           B[i]   = B[i]   + E*B[j];
       end
       println("");
       println("Passo ",j);
       println("--------------------------");
       PrintLinearSystem(A,B,PI,PD);
   end
   println("");
   println("---------------------------------------------------------------");
   println("P R O C E S S O   D E   S U B S T I T U I Ç Ã O   R E V E R S A");
   println("---------------------------------------------------------------");
   println("\nDado o Sistema Triangular Superior Ux=B\n");
   PrintLinearSystem(A,B,PI,PD);
   println("\n--------------------------");
   for i = n:-1:1
       S = 0;
       for j = i+1:n
           S = S + A[i,j]*X[j];
       end
       X[i] = (B[i] - S)/A[i,i];
       println("Passo ",n-i+1,": X[",i,"] = ",RealFormat(X[i],PI,PD));
       println("--------------------------");
   end
end
#----------------------------------------------------------------------
# Simulação da EliminacaoGaussiana baseado em Matriz Triangular Inferior
#----------------------------------------------------------------------
function Teste_EliminacaoGaussiana_TI(AA,BB,PI,PD)
   println("---------------------------------------------------------------------------");
   println("T R A N S F O R M A Ç Ã O   P O R   E L I M I N A Ç Ã O   G A U S S I A N A");
   println("---------------------------------------------------------------------------");
   println("\nDado o Sistema Linear Ax=B\n");
   A = copy(AA);
   B = copy(BB);
   PrintLinearSystem(A,B,PI,PD);
   m,n = size(A);
   X = zeros(n);
   for j = n:-1:2
       for i = 1:j-1
           E = -1*(A[i,j]/A[j,j]);
           A[i,:] = A[i,:] + E*A[j,:];
           B[i]   = B[i]   + E*B[j];
       end
       println("");
       println("Passo ",n-j+1);
       println("--------------------------");
       PrintLinearSystem(A,B,PI,PD);
   end
   println("");
   println("-------------------------------------------------------------");
   println("P R O C E S S O   D E   S U B S T I T U I Ç Ã O   D I R E T A");
   println("-------------------------------------------------------------");
   println("\nDado o Sistema Triangular Inferior Lx=B\n");
   PrintLinearSystem(A,B,PI,PD);
   println("\n--------------------------");
   for i = 1:n
       S = 0;
       for j = 1:i-1
           S = S + A[i,j]*X[j];
       end
       X[i] = (B[i] - S)/A[i,i];
       println("Passo ",i,": X[",i,"] = ",RealFormat(X[i],PI,PD));
       println("--------------------------");
   end
end
#-------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------
# Simulação do método Gauss-Jordan baseado em Matriz Diagonal
#----------------------------------------------------------------------

function Teste_GaussJordan(AA,BB,PI,PD)
   println("-----------------------------------------------------------");
   println("T R A N S F O R M A Ç Ã O   P O R   G A U S S - J O R D A N");
   println("-----------------------------------------------------------");
   println("\nDado o Sistema Linear Ax=B\n");
   A = copy(AA);
   B = copy(BB);
   PrintLinearSystem(A,B,PI,PD);
   m,n = size(A);
   X = zeros(n);
   for j = 1:n
       for i = 1:n
           if i!=j
              E = -1*(A[i,j]/A[j,j]);
              A[i,:] = A[i,:] + E*A[j,:];
              B[i]   = B[i]   + E*B[j];
           end
       end
       println("");
       println("Passo ",j);
       println("--------------------------");
       PrintLinearSystem(A,B,PI,PD);
   end
   println("");
   println("-------------------------------------------------------------");
   println("S U B S T I T U I Ç Ã O   D I R E T A");
   println("-------------------------------------------------------------");
   println("\nDado o Sistema Triangular Inferior Lx=B\n");
   PrintLinearSystem(A,B,PI,PD);
   println("\n--------------------------");
   for i = 1:n
       X[i] = B[i]/A[i,i];
       println("Passo ",i,": X[",i,"] = ",RealFormat(X[i],PI,PD));
       println("--------------------------");
   end
end
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
