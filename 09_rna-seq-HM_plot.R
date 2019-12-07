# Carregando a biblioteca xlsx para lidar com arquivos do Excel
library("xlsx")
# Carregando a biblioteca da função para gerar o HeatMap - heatmap.2()
library("gplots")
# Carregando a biblioteca com a função para criar a paleta de cores
library("RColorBrewer")

# Equivalente ao comando "pwd" no linux
getwd()
# Equivalente ao comando "cd" no linux
setwd("/state/partition1/lamoroso/")

getwd()

# Carrega arquivo texto separado por TAB ("\t") com cabeçalho e sem transformar 
# strings em fatores
deg.df <- read.delim(file="./output/13_trinity_results/abundance/DEG/SAMPLEA-SAMPLEB.txt", 
                     sep="\t", 
                     header=TRUE, 
                     stringsAsFactors = FALSE)

# exibir somente as 2 primeiras linhas para checagem
head(deg.df,2)

# Dimensão do data.frame, onde o primeiro valor refere-se à quantidade de linhas e
# o segundo valor à quantidade de colunas
dim(deg.df)

# selecionar um subconjunto de linhas do data.frame deg.df em um outro data.frame
# ("subset.deg.df") contendo somente os genes que possuírem logFC >= 2 ou logFC <= -2,
# além de possuírem um p-valor corrigido ("FDR") <= 0.05
subset.deg.df <- subset(deg.df, (((logFC >= 2) | (logFC <= -2)  ) & (FDR <=0.05))  )

dim(subset.deg.df)

head(subset.deg.df,2)

# criando um vetor de nomes de colunas selecionadas,
# ou seja, colunas que não são "X", "logFC", "PValue" ou "FDR",
# as quais, sabemos previamente que contém os valores de expressão.
# Os nomes das colunas serão ordenados depois de secionados com a 
# função "setdiff", a qual obtém
# um conjunto de valores (vetor) coma a diferença entre dois conjuntos de
# valores: o de todas as colunas ("colnames(subset.deg.df)") menos os nomes
# das colunas indesejadas.

sel_columns <- sort(
  setdiff( colnames(subset.deg.df), 
           c("X", "logFC", "PValue", "FDR") 
  )
)

# criando uma nova matriz "expression_data" que contém somente as colunas
# selecionadas
expression_data <- as.matrix( subset.deg.df[,sel_columns] )

# Atribuindo para os nomes das linhas da matriz o conteúdo da coluna "X",
# onde sabemos previamente que contém o identificador dos genes
rownames(expression_data) <- subset.deg.df$X

head(expression_data,2)

# Função para gravar o data.frame em um arquivo Excel, na planilha 
# de nome "DEGS_SAMPLEA-SAMPLEB"
write.xlsx(subset.deg.df,  
           file="./output/13_trinity_results/abundance/DEG/SAMPLEA-SAMPLEB.xlsx",
           sheetName="DEGS_SAMPLEA-SAMPLEB"
)


# cria uma paleta personalizada de 299 cores do vermelho ao verde,
# passando pelo amarelo
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

# define as quebras das cores manualmente para transição de cores
col_breaks = c(seq(-1,0,length=100),        # for red
               seq(0.01,0.8,length=100),    # for yellow
               seq(0.81,1,length=100))      # for green


# criação de uma imagem de tamanho 5 x 5 polegadas
png("./output/13_trinity_results/abundance/DEG/heatmap_DEGS_SAMPLEA-SAMPLEB.png",    # cria arquivo do tipo PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels de largura
    height = 5*300,       # 5 x 300 pixels de altura
    res = 300,            # 300 pixels por polegada
    pointsize = 8)        # tamanho da fonte

heatmap.2(expression_data, # a matriz com os valores de expressão
          main = "Correlation", # Título do HeatMap
          density.info="none",  # desabilita o gráfico de densidade dentro da legenda
          trace="none",         # desabilita as linhas dentro do HeatMap
          margins =c(12,12),     # definiação das margens no entorno do gráfico
          col=my_palette,       # nome do objeto contendo a paleta de cores criada anteriormente
          breaks=col_breaks,    # pontos de quebra para a transição de cores
          dendrogram="both",    # desenhar dendrograma para linhas e colunas
          distfun = function(x) as.dist(1-cor(t(x))), # distância baseada em correlação
          hclustfun = function(x) hclust(x, method="centroid") # método de ligação pelo centróide
)

dev.off()               # fecha o arquivo da imagem PNG
