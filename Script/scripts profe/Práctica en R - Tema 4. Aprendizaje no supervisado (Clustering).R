rm(list=ls())

path <- "/Users/vic/Library/CloudStorage/GoogleDrive-vdelaopascual@gmail.com/Mi unidad/MU en Bioinformática (UNIR 2023)/Actividades (AAAA-MM-DD)/Actividad 1_junio2024 update"

setwd(path)

df <- read.csv("Dataset expresión genes.csv")


library(dplyr)
df_genes <- df %>% dplyr::select(starts_with("AQ_"))
str(df_genes)


is.na(colSums(df_genes)) # ver si hay missing
df_genes_scale <- scale(df_genes)  # Normalización z-score



#### ---- Clustering no jerárquico con kmeans ----
library(factoextra)
library(stats)


kmeans.result <- kmeans(df_genes_scale, centers = 2, iter.max = 100, nstart = 25)
fviz_cluster(kmeans.result, df_genes_scale, xlab = '', ylab = '') +
  ggtitle("Cluster plot, centers = 2", subtitle = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, margin = margin(b = -10)))


kmeans.result <- kmeans(df_genes_scale, centers = 3, iter.max = 100, nstart = 25)
fviz_cluster(kmeans.result, df_genes_scale, xlab = '', ylab = '') +
  ggtitle("Cluster plot, centers = 3", subtitle = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, margin = margin(b = -10)))



# n optimo de clusters
fviz_nbclust(df_genes_scale, kmeans, method = "wss") +
  ggtitle("Optimal number of clusters", subtitle = "") +
  theme_classic()






#### ---- Clustering jerarquico aglomerativo (de abajo a arriba): datos con poca cantidad de observaciones ----
library(ggdendro)
library(cluster)

# Calcular la matriz de distancia
dist_matrix <- dist(df_genes_scale)

# Se ejecuta el algoritmo de clusterización jerárquica aglomerativa
hclust_model_single <- hclust(dist_matrix, method = "single") # agrupa los clusters usando la distancia entre los puntos más CERCANOS
hclust_model_complete <- hclust(dist_matrix, method = "complete") # agrupa los clusters usando la distancia entre los puntos más ALEJADOS
hclust_model_average <- hclust(dist_matrix, method = "average") # agrupa los clusters usando el PROMEDIO de todas las distancias entre los puntos de ambos clusters
hclust_model_ward <- hclust(dist_matrix, method = "ward.D") # agrupa los clusters tratando de que sean lo más COMPACTOS posible minimizando la dispersión interna


library(ggplot2) 
library(factoextra)
library(cluster)

colors <- rainbow(5)

# single: conecta puntos cercanos y puede encontrar clusters largos y delgados, PERO a veces conecta muchos puntos formando cadenas, lo que puede ser poco útil
clust_single <- fviz_dend(hclust_model_single, 
                          cex = 0.5,
                          k = 5,
                          palette = colors,
                          main = "Single",
                          xlab = "Índice de Observaciones",
                          ylab = "Distancia") + 
  theme_classic()

# complete: crea clusters compactos y bien definidos, PERO es sensible a puntos extremos (outliers), que pueden distorsionar los clusters
clust_complete <- fviz_dend(hclust_model_complete, 
                          cex = 0.5,
                          k = 5,
                          palette = colors,
                          main = "Complete",
                          xlab = "Índice de Observaciones",
                          ylab = "Distancia") + 
  theme_classic()

# complete: encuentra un equilibrio entre single y complete, PERO a veces no es tan bueno para clusters con tamaños muy diferentes
clust_average <- fviz_dend(hclust_model_average, 
                          cex = 0.5,
                          k = 5,
                          palette = colors,
                          main = "Average",
                          xlab = "Índice de Observaciones",
                          ylab = "Distancia") + 
  theme_classic()

# Ward: crea clusters redondeados y homogéneos similares a los que genera k-means, PERO no funciona tan bien si los clusters tienen formas raras o tamaños muy diferentes
clust_ward <- fviz_dend(hclust_model_ward, 
                          cex = 0.5,
                          k = 5,
                          palette = colors,
                          main = "Ward",
                          xlab = "Índice de Observaciones",
                          ylab = "Distancia") + 
  theme_classic()

library(gridExtra)
grid.arrange(clust_single, clust_complete, clust_average, clust_ward, nrow = 2)


# single: datos con patrones lineales -> "Creo que los datos tienen estructuras locales fuertes, y estoy más interesado en cómo se conectan los puntos cercanos."
# complete: datos donde esperas clusters compactos y bien definidos -> "Los datos están agrupados en regiones claramente separadas y no quiero que un punto extremo distorsione el análisis"
# average: cuando no tienes una forma clara en mente para los clusters -> "Espero una mezcla de clusters compactos y algo más dispersos, pero quiero un balance entre lo local y lo global."
# ward.D: datos en los que esperas clusters compactos y homogéneos -> "Los datos deberían agruparse en clusters compactos con baja variabilidad interna."



df$cluster_single <- as.factor(cutree(hclust_model_single, k = 5))
df$cluster_complete <- as.factor(cutree(hclust_model_complete, k = 5))
df$cluster_average <- as.factor(cutree(hclust_model_average, k = 5))
df$cluster_ward <- as.factor(cutree(hclust_model_ward, k = 4))





#### ---- Clustering jerarquico divisivo (de arriba a abajo): datos con muchas observaciones ----
library(ggplot2) 
library(factoextra)
library(cluster)

# Implementación del clustering divisivo
diana_euclidean <- diana(df_genes_scale, metric = "euclidean", stand = FALSE) # ideal para datos donde las distancias más pequeñas
diana_manhattan <- diana(df_genes_scale, metric = "manhattan", stand = FALSE) # ideal para datos con diferentes escalas o datos categóricos


colors <- rainbow(5)
clust_diana_euclidean <- fviz_dend(diana_euclidean, 
                                   cex = 0.5, 
                                   k = 5,
                                   palette = colors, 
                                   main = 'Euclidean',
                                   xlab = "Índice de Observaciones",
                                   ylab = "Distancia") + 
  theme_classic()


colors <- rainbow(5)
clust_diana_manhattan <- fviz_dend(diana_manhattan, 
                                   cex = 0.5, 
                                   k = 5,
                                   palette = colors, 
                                   main = 'Manhattan',
                                   xlab = "Índice de Observaciones",
                                   ylab = "Distancia") + 
  theme_classic()


grid.arrange(clust_diana_euclidean, clust_diana_manhattan, nrow = 2)









rm(list=ls())

# Cargar librerías
library(ggplot2)
library(cluster)

# Ejemplo de coordenadas (representación de plantas, clientes, tiendas, etc.)
data <- data.frame(
  x = c(2, 4, 6, 10),  # Coordenadas x (pueden ser precios, calificaciones, etc.)
  y = c(3, 4, 7, 10)   # Coordenadas y
)

# Calcular la matriz de distancias (usando Euclidiana)
dist_matrix <- dist(data)
dist_matrix

# Realizar el clustering con diferentes métodos de enlace
hclust_single <- hclust(dist_matrix, method = "single")
hclust_complete <- hclust(dist_matrix, method = "complete")
hclust_average <- hclust(dist_matrix, method = "average")
hclust_ward <- hclust(dist_matrix, method = "ward.D2")

hclust_single$order
hclust_complete$order
hclust_average$order
hclust_ward$order

# Graficar resultados
par(mfrow = c(2, 2))  # Dividir la ventana en 2x2 para graficar múltiples
plot(hclust_single, main = "Single Linkage")
plot(hclust_complete, main = "Complete Linkage")
plot(hclust_average, main = "Average Linkage")
plot(hclust_ward, main = "Ward's Linkage")




