# Trabalho de álgebra linear computacional
## Professor João Paixão
### Programa de pós graduação em informática (PPGI) 
#### Universidade Federal do Rio de Janeiro (UFRJ)


##### A ideia deste trabalho é criar uma biblioteca do 0 que implemente na linguagem Julia as funções aprendidas na disciplina sobre álgebra linear.



**Não sabe ainda como instalar o Jupyter notebook? :sweat_smile:**

[Segue esse tutorial aqui que não tem erro](https://medium.com/horadecodar/como-instalar-o-jupyter-notebook-windows-e-linux-20701fc583c)


**Vixi, mas não tem Julia no Jupyter! :scream:**
O Jupyter sozinho não vem com Julia, precisa de uma extensão.

[Segue esse tutorial aqui que não tem erro](https://datatofish.com/add-julia-to-jupyter/)

 
Não instala o anaconda, ele é um misto de vários programas pra quem trabalha com análise de dados e python. Para nós, vai ser uma bazuca pra matar mosquito. :honeybee: :gun:

**Léo, que raios é git e github?**
Separei um artigo super legal que explica de uma maneira extremamente simples o que é e como configurar e começar a utilizar. Quaisquer dúvidas, entrem em contato! 

[Segue esse tutorial aqui que não tem erro](https://medium.com/reprogramabr/git-e-github-por-onde-come%C3%A7ar-ca88a783c223)


**Como vai funcionar?**

Clone este repositório com o comando `git clone`, insira as suas funções e realize um `pull request`. 

Qualquer dúvida, só perguntar! :blush:

**Lista de funcionalidades que a biblioteca deve ter**

| Especificações                                               |                                           Entrada | Saída                            | Responsável | Status |
|--------------------------------------------------------------|--------------------------------------------------:|----------------------------------|-------------|--------|
| Produto interno                                              | 2 vetores u, v                                    | escalar c                        | leo         |	:white_check_mark:	 | 
| Posto da matriz                                              | 1 matriz                                          | inteiro                          | nomair      |	:x:	 |
| Comprimento de um vetor (norma do vetor)                     | 1 vetor                                           | escalar                          | leo         |	:x:	 |
| Ângulo entre dois vetores                                    | 2 vetores u, v                                    | ângulo theta (escalar)           | Nomair      |	:x:	 |
| Se 2 vetores são perpendiculares                             | 2 vetores: u,v                                    | booleano                         | Nomair      |	:x:	 |
| Se vetores são ortogonais                                    | vetores                                           | booleano                         | leo         |	:x:	 |
| Se duas retas são ortogonais                                 | 2 retas: 2 vetores: u,v                           | booleano                         | waldir      |	:x:	 |
| Se uma vetor é ortogonal a um plano                          | 1 reta, 1 plano                                   | booleano                         | leo         |	:x:	 |
| Eliminação Gaussiana                                         | 1 matriz quadrada A e 1 vetor coluna b            | Um vetor x                       | gastao      |	:x:	 |
| Verificar se uma matriz é triangular                         | 1 matriz                                          | booleano                         | gabi        |	:x:	 |
| Verificar se uma matriz é ortogonal                          | 1 matriz A                                        | booleano                         | gabi        |	:x:	 |
| Verficar se um conjunto de vetores forma uma base ortonormal | Conjunto de vetores                               |                                  | Gabi        |	:x:	 |
| Projeção de um vetor em uma reta                             | Vetor u e vetor direção da reta v                 | Vetor projetado a                | waldir      |	:x:	 |
| Projeção de um vetor em um plano                             | 1 vetor, 1 plano                                  |                                  | waldir      |	:x:	 |
| Verificar se uma matriz é inversível                         | 1 matriz                                          | booleano                         |             |	:x:	 |
| Gauss-jordan                                                 | 1 matriz e 1 vetor                                | 1 vetor x                        | pablo       |	:x:	 |
| Decomposiçao LU                                              | Uma matriz A                                      | Duas matrizes triangulares L e U | pablo       |	:x:	 |
| Substituiçao reversa (ou recursiva)                          | 1 matriz triangular superior A e 1 vetor coluna b | vetor x                          | Gabi        |	:x:	 |
| Substituição normal (ou recursiva)                           | 1 matriz triangular inferior A e 1 vetor coluna b | vetor x                          | gastao      |	:x:	 |
| Mínimos quadrados                                            | Matriz A e vetor b                                | Vetor x                          |             |	:x:	 |
| Decomposição QR                                              | 1 matriz                                          | 2 matrizes                       | Gastão      |	:x:	 |
| SVD (só o primeiro vetor)                                    | Conjunto de pontos: matriz                        | vetor v                          | gastao      |	:x:	 |
| Determinante                                                 | Uma matriz quadrada                               | escalar                          | pablo       |	:x:	 |
| Produto vetorial                                             | 2 vetores p, q                                    | vetor r                          |             |	:x:	 |
| Pseudo inversa                                               | matriz A                                          | matriz A*                        | gabi        |	:x:	 |
| Achar a menor solução de um sistema linear                   | matriz A, vetor b                                 | vetor x*                         |             |	:x:	 |
| Achar infinitas soluções de um sistema linear                | matriz A, vetor b                                 | ??                               |             |	:x:	 |
| Achar nucleo soluçao sistema linear                          | 1 matriz                                          | conjunto vetores                 |             |	:x:	 |
| Achar os b que o sistema Ax=b tem solução                    | lista de variáveis e relações entre elas          | vetor b                          |             |	:x:	 |
| Sistema dinâmicos lineares                                   | Matriz A, um vetor x e um inteiro k               | Uma matriz C                     | gastao      |	:x:	 |
| Modelagem                                                    | lista de variáveis e relações entre elas          | uma matriz de solução            |             |	:x:	 |
| Normas                                                       | 1 vetor b                                         | escalar c                        | Waldir      |	:x:	 |
| Produto de matrizes (ou recursiva)                           | Duas matrizes A e B                               | Uma matriz C                     | Waldir      |	:x:	 |
| Soma de matrizes (ou recursiva)                              | Duas matrizes A e B                               | Uma matriz X                     | Waldir      |	:x:	 |
| Encontrar Inversa                                            |                                                   |                                  |             |	:x:	 |
| Testes: Matrizes do problemas de modelagem, aleatórias       | ??                                                |                                  |             |	:x:	 |
| Gauss-Seidel                                                 | 1 matriz A, 1 vetor b                             | Um vetor x                       | pablo       |	:x:	 |
| Gauss-Jacobi                                                 | 1 matriz A, 1 vetor b                             | Um vetor x                       | pablo       |	:x:	 |
| Visualizações, input-output                                  |                                                   |                                  | pablo       |	:x:	 |
| Transposta                                                   |                                                   |                                  |             |	:x:	 |