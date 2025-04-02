## ----include=FALSE------------------------------------------------------------
# Feature data, el cual será nuestra matriz
features <- read.csv("C:/Users/verag/OneDrive/Escritorio/UOC/Análisis de datos ómicos/Garmendia-Dulfova-Vera-PEC1/dataset/features.csv", sep = ";")
cat("Estructura de features:\n")
str(features)
# Sample data
metadata <- read.csv("C:/Users/verag/OneDrive/Escritorio/UOC/Análisis de datos ómicos/Garmendia-Dulfova-Vera-PEC1/dataset/metadata.csv", sep = ";")
cat("\nEstructura de metadata:\n")
str(metadata)
# Metabolite data 
metabolites <- read.csv("C:/Users/verag/OneDrive/Escritorio/UOC/Análisis de datos ómicos/Garmendia-Dulfova-Vera-PEC1/dataset/metaboliteNames.csv", sep = ";")
cat("\nEstructura de metabolites:\n")
str(metabolites)


## ----include=FALSE------------------------------------------------------------
# Convertimos el conjunto de datos features en una matriz
features_matrix <- as.matrix(features)
# y nos aseguramos de que dicha conversión se ha realizado correctamente.
str(features_matrix)


## -----------------------------------------------------------------------------
dim(features_matrix) # matriz
dim(metadata) # sample metadata
dim(metabolites) # feature metadata


## ----include=FALSE------------------------------------------------------------
# Convertimos la columna 'PubChem' en rownames
rownames(metabolites) <- metabolites$PubChem
# Eliminamos la columna 'PubChem' 
metabolites$PubChem <- NULL

# Convertimos la columna 'ID' en rownames
rownames(metadata) <- metadata$ID
# Eliminamos la columna 'ID'
metadata$ID <- NULL

# Ordenamos metabolites por PubChem
metabolites_ordered <- metabolites[order(rownames(metabolites)), ]

# Ordenamos features_matrix por los nombres de las filas
features_matrix_ordered <- features_matrix[order(rownames(features_matrix)), ]


## ----include=FALSE------------------------------------------------------------
library(SummarizedExperiment)


## -----------------------------------------------------------------------------
stopifnot(rownames(features_matrix_ordered) == metabolites_ordered$PubChem)
stopifnot(colnames(features_matrix_ordered) == metadata$ID)

# Creamos el objeto SummarizedExperiment
se <- SummarizedExperiment(assays = list(counts = features_matrix_ordered),
                           colData = metadata,
                           rowData = metabolites_ordered)
se


## ----include=FALSE------------------------------------------------------------
# Guardamos el objeto "se" en formato .Rda
save(se, file = "C:/Users/verag/OneDrive/Escritorio/UOC/Análisis de datos ómicos/Garmendia-Dulfova-Vera-PEC1/ObjetoSummarizedExperiment.Rda")


## -----------------------------------------------------------------------------
# Accedemos a la matriz del objseto se y miramos si hay valores NA
sum(is.na(assay(se)))
# Como hay muchos valores NA, los imputamos con la media de cada columna
assay(se) <- apply(assay(se), 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
# Verificamos que no queden valores NA
sum(is.na(assay(se)))


## -----------------------------------------------------------------------------
# Asignamos colores a cada tratamiento
groupColors <- c("Baseline" = "blue", "Apple" = "green", "Cranberry" = "red")

# Creamos un vector de colores para cada muestra según su tratamiento
col <- groupColors[colData(se)$Treatment]


## -----------------------------------------------------------------------------
# Filtramos las muestras por tratamiento
baseline <- se[, se$Treatment == "Baseline"]
apple <- se[, se$Treatment == "Apple"]
cranberry <- se[, se$Treatment == "Cranberry"]

cat("Baseline\n")
round(apply(assay(baseline[1:5,]),1, summary))
cat("\nApple\n")
round(apply(assay(apple[1:5,]),1, summary))
cat("\nCranberry\n")
round(apply(assay(cranberry[1:5,]),1, summary))


## -----------------------------------------------------------------------------
boxplot(assay(se), col = col, main="Distribucion de los valores de expresión")


## -----------------------------------------------------------------------------
boxplot(log(assay(se)+1), col = col, main="Distribución de los valores de metabolitos (transformación logarítmica)")


## -----------------------------------------------------------------------------
# Realizamos el Análisis de Componentes Principales
pca <- prcomp(t(log(assay(se)+1)), scale. = TRUE)


## -----------------------------------------------------------------------------
# Calculamos la varianza explicada por cada componente principal
loads <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 1)


## -----------------------------------------------------------------------------
# Creamos el gráfico de dispersión de los dos primeros componentes principales
xlab <- c(paste("PC1", loads[1], "%"))
ylab <- c(paste("PC2", loads[2], "%"))
plot(pca$x[, 1:2], xlab = xlab, ylab = ylab, col = col, 
     main = "Principal components (PCA)")

names2plot <- rownames(colData(se))

text(pca$x[,1],pca$x[,2], names2plot, pos=3, cex=.6)


## -----------------------------------------------------------------------------
# Definimos el método de cálculo de distancia
distmeth <- c("euclidean")
Distan <- dist(t(assay(se)), method=distmeth) # Matriz de distancia entre muestras
# Definimos el método de agrupamiento
treemeth <- c("average")
hc <- hclust(Distan, method=treemeth) 
plot(hc)

