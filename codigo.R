# HOLA MUNDO

saludo <- "¡Hola Mundo!"

saludo

# USANDO FUNCIONES head()

df <- iris

head(df)

# USANDO FUNCIONES pairs()

df <- iris

pairs(df, col = df$Species)


# if, else

variable_1 <- 5

if (variable_1 >10) {
  print('la variable_1 es mayor que 10')
} else {
  print('la variable_1 NO es mayor que 10 ')
}


# for
estaciones <- c('primavera', 'verano', 'otoño', 'invierno')

for(estacion in estaciones){
  print(toupper(estacion))
}

# VARIABLES SIMPLES

numero <- 1.2
numero

entero <- 10
entero

texto <- 'esto es un texto'
texto

complejo <- 3 + 2i
complejo

logical <- 5 == 10
logical

# VECTORES

vector_numeros <- c(2, 4, 6, 8)
vector_numeros

seq_numeros <- 10:20
seq_numeros

comb_vect <- c(vector_numeros, seq_numeros)
comb_vect

# Aritmetic de Vectores

# Un vector
vector_numeros <- seq(10, 50, by = 5)
vector_numeros

vector_numeros + 1

vector_numeros / 2

#Dos vectores

vector_numeros_1 <- c(1, 2, 3, 4)
vector_numeros_2 <- c(1, 1, 0, 0)

vector_numeros_1 - vector_numeros_2

vector_numeros_1 * vector_numeros_2

#Vector de numeros y vector de texto

letras <- c('A', 'B', 'C', 'D')
numeros <- c(1, 2, 3, 4)

paste0(letras, numeros)

paste0(numeros, letras)

# Elementos de un vector

numeros <- 10:20
numeros

numeros[4]
numeros[-4]
numeros[c(4:8)]
numeros[-c(4:8)]
numeros[c(2, 10)]

numeros[1] <- 0
numeros

numeros[c(2, 10)] <- 200
numeros


numeros <- 10:20
numeros

numeros <- numeros[-c(4:8)]
numeros

numeros <- c(numeros, c(1, 2, 3))
numeros


#Mascara logica

numeros <- 10:20
numeros

numeros < 15
numeros[numeros < 15]

mascara_logica <- numeros == 10

mascara_logica

!mascara_logica

numeros[!mascara_logica]



#vector indices
numeros <- runif(10, min = 1, max = 1000)
numeros
numeros[c(1, 5, 10)]

#mascara_logica
numeros[numeros < 500]



# DATA FRAME

df <- iris

dim(df)

head(df)

summary(df)


df[ , ]

## Structure of an Arbitrary R Object

str(df)

df$Sepal.Length

df[, c('Sepal.Length')]
df[, c(1)]

### Seleccionar filas y columnas

df[5:10, c('Sepal.Length', 'Species')]

# Mascara Logica

head(df)

df[df$Species == 'versicolor' & df$Sepal.Length > 6.5, ]


# Add column

df$Ratio <- df$Sepal.Length / df$Sepal.Width
head(df)

df[df$Species == 'setosa', c('Ratio')] <- 0
df

# DEL Column
df$Ratio <- NULL
head(df)


# Crear un data frame

df.ratio <- data.frame(
  Species = df$Species,
  Sepal.Ratio = df$Sepal.Length / df$Sepal.Width,
  Petal.Ratio = df$Petal.Length / df$Petal.Width
)
head(df.ratio)

#########################################################################
#########################################################################
##### ggplot2 ####

# install.packages('ggplot2')

library(ggplot2)
df <- diamonds
summary(df)

###### Scatterplot

ggplot(df, aes(x = carat, y = price)) +
  geom_point()

# añadir color y temas
ggplot(df, aes(x = carat, y = price)) +
  geom_point(color = 'steelblue') +
  theme_bw()

# colorear cada punto por valor de cut
ggplot(df, aes(x = carat, y = price)) +
  geom_point(aes(color = cut)) +
  theme_bw() +
  theme(legend.position = 'top')

# facet_wrap
ggplot(df, aes(x = carat, y = price)) +
  geom_point(aes(color = cut)) +
  theme_bw() +
  theme(legend.position = 'top') +
  facet_wrap(~ cut)

####### Boxplot
ggplot(df, aes(y = price)) +
  geom_boxplot()

# añadir cut 
ggplot(df, aes(x = cut, y = price)) +
  geom_boxplot(color = 'navy') +
  theme_bw()

# fondo cut
ggplot(df, aes(x = cut, y = price)) +
  geom_boxplot(aes(fill = cut), color = 'navy') +
  theme_bw() +
  theme(legend.position = 'top')

# Añadir variable clarity
ggplot(df, aes(x = cut, y = price)) +
  geom_boxplot(aes(fill = clarity), color = 'navy') +
  theme_bw() +
  theme(legend.position = 'top')

# facet_wrap
ggplot(df, aes(x = cut, y = price)) +
  geom_boxplot(aes(fill = clarity), color = 'navy') +
  theme_bw() +
  theme(legend.position = 'top') +
  facet_wrap(~ clarity)

# Bar chart
ggplot(df, aes(x = clarity)) +
  geom_bar(color = 'navy', fill = 'white') +
  theme_bw()

# añadir color clarity
ggplot(df, aes(x = clarity)) +
  geom_bar(aes(color = clarity, fill = clarity)) +
  theme_bw() +
  theme(legend.position = 'top')

#añadir color cut
ggplot(df, aes(x = clarity)) +
  geom_bar(aes(color = cut, fill = cut)) +
  theme_bw() +
  theme(legend.position = 'top')

# position dodge
ggplot(df, aes(x = clarity)) +
  geom_bar(aes(color = cut, fill = cut), position = position_dodge()) +
  theme_bw() +
  theme(legend.position = 'top')


#####################################################################3
#####################################################################
##### RNA-SEQ

df <- read.csv('data/mcf_count_table.csv')

dim(df)


metadata <- read.csv('data/mcf_metadata.csv')
dim(metadata)


## Fenotipos seleccionados

fenotipos <- c('L1Larvae', 'L2Larvae')

#### df.fenotipos <- metadata[metadata$stage == 'L1Larvae' | 
####                           metadata$stage == 'L2Larvae'
####                            , c('sample.id', 'stage')]

df.fenotipos <- metadata[metadata$stage %in% fenotipos, c('sample.id', 'stage')]

df.fenotipos$stage <- as.factor(df.fenotipos$stage)

df.fenotipos


#### seleccionar las muestras de los fenotipos seleccionados
names(df)

df.sub <- df[, df.fenotipos$sample.id]

head(df.sub)

rownames(df.sub) <- df$gene

head(df.sub)


## ELIMINAR genes con pocas counts

head(rowSums(df.sub), 30)

keep <- rowSums(df.sub) >= 10
head(keep, 30)


df.sub  <- df.sub [keep, ]
head(df.sub)


#### DESeq2
#### install.packages("BiocManager")
#### BiocManager::install("DESeq2")


library(DESeq2)
coldata <- data.frame(stage = df.fenotipos$stage)
rownames(coldata) <- df.fenotipos$sample.id
head(coldata)


## Objeto DESeqDataSet
dds <- DESeqDataSetFromMatrix(
                      countData = df.sub,
                      colData = coldata,
                      design = ~stage
              )

dds


dds <- DESeq(dds)
dds


res <- results(dds, contrast = c('stage', 'L1Larvae', 'L2Larvae'))
res


### Guardar resultados en CSV

write.csv(res, 'data/DGE.csv')

rld <- assay(rlog(dds, blind = FALSE))
head(rld)

write.csv(rld, 'data/normalized_counts.csv')



##################################################
######## Volcano PLot

dge <- read.csv('data/DGE.csv')

#### ELIMINAR LOS GENES CON NA values
dge <- dge[complete.cases(dge), ]

head(dge)


### Valores de corte
padj.cutoff <- 0.0000001
lfc.cutoff <- log2(2)


### clasificar genes
dge$class <- 'none'

dge[dge$log2FoldChange >= lfc.cutoff & dge$padj <= padj.cutoff, c('class')] <- 'UP'
dge[dge$log2FoldChange <= -1 * lfc.cutoff & dge$padj <= padj.cutoff, c('class')] <- 'DOWN'

head(dge)


ggplot(dge, aes(x = log2FoldChange, y = -1 * log10(padj))) +
  geom_point() +
  theme_bw()

### añadir rectar de valores de corte
ggplot(dge, aes(x = log2FoldChange, y = -1 * log10(padj))) +
  geom_point() +
  geom_hline(yintercept = -1 * log10(padj.cutoff ), linetype="dashed", 
             color = "firebrick", linewidth = 0.5) +
  geom_vline(xintercept=c(-1 * lfc.cutoff ,lfc.cutoff ), linetype="dashed", 
             color = "firebrick", linewidth = 0.5) +
  theme_bw()


ggplot(dge, aes(x = log2FoldChange, y = -1 * log10(padj))) +
  geom_point(aes(color = class)) +
  geom_hline(yintercept = -1 * log10(padj.cutoff ), linetype="dashed", 
             color = "firebrick", linewidth = 0.5) +
  geom_vline(xintercept=c(-1 * lfc.cutoff ,lfc.cutoff ), linetype="dashed", 
             color = "firebrick", linewidth = 0.5) +
  theme_bw()



### TOP GENES

dge <- dge[order(dge$log2FoldChange, decreasing = TRUE), ]
genes.up <- dge[1:40, c('X')]
genes.down <- dge[(nrow(dge) - 30):nrow(dge), c('X')]

genes.up
genes.down

df.norm <- read.csv('data/normalized_counts.csv')
head(df.norm)


df.norm <- df.norm[df.norm$X %in% c(genes.up, genes.down), ]


###### HEATMAP
#install.packages('pheatmap')
library(pheatmap)

pheatmap(df.norm[, 2:10])

pheatmap(df.norm[, 2:10], 
         labels_row = df.norm$X
         )

pheatmap(df.norm[, 2:10], 
         labels_row = df.norm$X,
         scale = 'row'
        )


pheatmap(df.norm[, 2:10], 
         labels_row = df.norm$X,
         scale = 'row',
         cutree_cols =  2,
         cutree_rows = 2,
         color = colorRampPalette(c("steelblue", "white", "firebrick3"))(100)
        )

pheatmap(df.norm[, 2:10], 
         labels_row = df.norm$X,
         scale = 'row',
         cutree_cols =  2,
         cutree_rows = 2,
         color = colorRampPalette(c("steelblue", "white", "firebrick3"))(100),
         annotation_col = coldata
)



##### ENSEMBL ID
### BiocManager::install("org.Dm.eg.db")

library(org.Dm.eg.db)

keytypes(org.Dm.eg.db)


genes.names <- mapIds(org.Dm.eg.db, 
                      keys = df.norm$X, 
                      keytype = 'ENSEMBL', 
                      column = 'SYMBOL')

genes.names

genes.names <- as.data.frame(genes.names)
genes.names

genes.names$X <- rownames(genes.names)
genes.names

df.norm.names <- merge(df.norm, genes.names, by = 'X')


pheatmap(df.norm.names[, 2:10], 
         labels_row = df.norm.names$genes.names,
         scale = 'row',
         cutree_cols =  2,
         cutree_rows = 2,
         color = colorRampPalette(c("steelblue", "white", "firebrick3"))(100),
         annotation_col = coldata
)
