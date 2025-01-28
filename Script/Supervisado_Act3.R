#las "-------------" separan los bloques de codigo para el rmd

rm(list=ls())

setwd("/Users/carme/OneDrive/Escritorio/R/Master_bioinformatica/Algoritmos/actividades/Act3_IA")

#Preparación del entorno de trabajo

#Cargar librerias
library(glmnet) # ElasticNet
library(tidyverse)
library(caret) # ML
library(rpart) # DT
library(rpart.plot) # DT plot
library(rattle) # DT plot
library(pROC) # ROC
library(PRROC) # PR-Curve
library(gridExtra) # juntar los gráficos
library(dplyr)
library(MASS)
library(class)
library(klaR)
library(rpart)

#Cargar dataset
classes <- read.csv("Datos/classes.csv", header=F, sep = ";")
gene_expression <- read.csv("Datos/gene_expression.csv", header=F, sep = ";", col.names = read_lines("Datos/column_names.txt"))


------------------------------------------------------------------------------------------
# Nombre de las columnas
  
  colnames(gene_expression) <- read_lines("Datos/column_names.txt")

# Asignar nombres a las columnas de classes

colnames(classes) <- c("Muestra", "Clase")

# Datos NA 
sumas <- colSums(gene_expression) 
columnascero <- names(sumas[sumas==0])
gene_expression <- gene_expression[, !names(gene_expression) %in% columnascero] 

# Unir archivos

dataframe <-cbind(gene_expression, classes)


# Mover columnas clase y muestra
dataframe <- dataframe %>%
  dplyr::select(Muestra, Clase, everything())

#Escalado de datos
data_scaled <- dataframe %>%
  mutate(across(where(is.numeric), scale))


------------------------------------------------------------------------------------------



# Dividir el conjunto de datos: entrenamiento (80%) y prueba (20%)
data_scaled$Clase <- as.factor(data_scaled$Clase)
levels(data_scaled$Clase)

set.seed(1995)

# Crear índices de partición
train_index <- createDataPartition(data_scaled$Clase, p = 0.8, list = FALSE)

# Crear conjuntos de entrenamiento y prueba
train_data <- data_scaled[train_index, ]
test_data <- data_scaled[-train_index, ]



#A partir de aqui son cosas para la division de los datos que me ha dicho chatgpt con otras cosas que he visto del profesor pero no lo usé para la actividad 2

# Asegurarse de que la variable 'Muestra' es un factor
train_data$Muestra <- factor(train_data$Muestra)
test_data$Muestra <- factor(test_data$Muestra, levels = levels(train_data$Muestra))

# Seleccionar las columnas numéricas para escalarlas
numerical_columns_train <- train_data[, sapply(train_data, is.numeric)]
scaled_train_data <- scale(numerical_columns_train)

# Convertir los datos escalados a un data.frame
scaled_train_data <- as.data.frame(scaled_train_data)

# Unir los datos escalados con la variable Clase y Muestra en el conjunto de entrenamiento
train_data <- cbind(scaled_train_data, Clase = train_data$Clase, Muestra = train_data$Muestra)

# Asegurarse de que el conjunto de datos de entrenamiento es un data.frame
train_data <- as.data.frame(train_data)

# Escalar las columnas numéricas del conjunto de prueba
numerical_columns_test <- test_data[, sapply(test_data, is.numeric)]
scaled_test_data <- scale(numerical_columns_test)

# Convertir los datos escalados a un data.frame
scaled_test_data <- as.data.frame(scaled_test_data)

# Unir los datos escalados con la variable Clase y Muestra en el conjunto de prueba
test_data <- cbind(scaled_test_data, Clase = test_data$Clase, Muestra = test_data$Muestra)

# Asegurarse de que el conjunto de datos de prueba es un data.frame
test_data <- as.data.frame(test_data)

# Convertir 'Clase' a factor en ambos conjuntos de datos
train_data$Clase <- as.factor(train_data$Clase)
test_data$Clase <- as.factor(test_data$Clase)

# Confirmar la estructura de los datos
str(train_data)
str(test_data)

------------------------------------------------------------------------------------------

###########################################
#Modelo LDA (Linear Discriminant Analysis)#
###########################################
library(MASS)

#Entrenar el modelo LDA
lda_model <- MASS::lda(Clase ~ ., data = train_data)

#Observar plano
lda_model$scaling

#Comentar segun salga: ejemplO:
#Si un coeficiente es positivo, significa que un aumento en esa variable aumentará el valor de la función discriminante.
#Si un coeficiente es negativo, significa que un aumento en esa variable disminuirá el valor de la función discriminante.

#perimeter1 (6.530946007) tiene un coeficiente muy grande en LD1, lo que significa que es una variable muy importante para separar las clases en la primera función discriminante.
#radius1 (-7.574833219)tiene un coeficiente negativo significativo en LD1, lo que indica que un aumento en radius1 tiende a desplazar la clasificación.


lda_pred <- predict(lda_model, newdata = train_data)
lda_pred$x

#Realizar predicciones sobre el conjunto de prueba
lda_predictions <- predict(lda_model, newdata = test_data)
lda_predictions$x

#Obtener la predicción (predicciones de la clase)
predicted_classes <- lda_predictions$class
predicted_classes
length(predicted_classes)

#Obtener las verdaderas etiquetas (las clases reales en el conjunto de prueba)
true_classes <- as.factor(test_data$Clase)
true_classes
length(true_classes)

#Crear la matriz de confusión (predicho test vs. real test
lda_conf_matrix <- confusionMatrix(predicted_classes, true_classes)
print(lda_conf_matrix)

#Comentar tabla segun salga: ejemplO:
#Acuracy: El modelo tiene una precisión del 94,69%. Esto significa que el 94,69% de las predicciones fueron correctas (el modelo acertó en la clasificación).
#Sensitivity: La sensibilidad mide cuántas veces el modelo predijo correctamente la clase B entre las veces que realmente era  clase B. En este caso, es 1.0000, lo que significa que el modelo detectó todas las instancias de  clase B correctamente (100%).
#Specificity: La especificidad mide cuántas veces el modelo predijo correctamente que no era  clase B. En este caso, 85.71% de las veces que no era  clase B, el modelo lo predijo correctamente.
#Pos Pred Value: Este valor indica el porcentaje de veces que el modelo predijo  clase B correctamente cuando hizo esa predicción. Aquí es 92,21%.
#Neg Pred Value: Mide cuántas veces el modelo predijo correctamente que no era  clase B cuando realmente no era  clase B. En este caso, el modelo predijo correctamente 100% de las veces que no era  clase B.
#Balanced Accuracy: es un promedio de la sensibilidad y especificidad. En este caso, es 92,86%, lo que sugiere que el modelo tiene un buen desempeño en la  clase B.



#Gráfico
predict(lda_model)

library(ggplot2)
lda.data <- cbind(train_data, predict(lda_model)$x)
ggplot(lda.data, aes(x = LD1, fill = Clase)) + 
  geom_density(alpha = 0.5) + 
  theme_classic()

#Comentar plot segun salga......... Hay una buena separación entre las distribuciones de ambas clases, aunque existe una pequeña zona de solapamiento cerca de LD1 = 0. Esto sugiere que el modelo LDA tiene un buen desempeño para diferenciar entre las dos clases, pero puede haber ligeros errores en la clasificación en esa región de solapamiento.


------------------------------------------------------------------------------------------

##############################################
#Modelo QDA (No linear Discriminant Analysis)#
##############################################
library(MASS)

#Entrenar modelo QDA
qda_model <- MASS::qda(Clase ~ ., data = train_data)


#Realizar predicciones sobre el conjunto de prueba
qda_predictions <- predict(qda_model, newdata = test_data)

#Obtener la predicción (predicciones de la clase)
predicted_classes <- qda_predictions$class
predicted_classes
length(predicted_classes)

#Obtener las verdaderas etiquetas (las clases reales en el conjunto de prueba)
true_classes <- as.factor(test_data$Clase)
true_classes
length(true_classes)

#Crear la matriz de confusión (trat predicho testing vs. trat real testing)
qda_conf_matrix <- confusionMatrix(predicted_classes, true_classes)
print(qda_conf_matrix)

#Comentar tabla segun salga: ejemplO:
#Acuracy: el modelo tiene una precisión del 97.35%. Esto significa que el 97.35%. de las predicciones fueron correctas (el modelo acertó en la clasificación).
#Sensitivity: la sensibilidad mide cuántas veces el modelo predijo correctamente la clase B entre las veces que realmente era clase B. En este caso es 97,18%.
#Specificity: la especificidad mide cuántas veces el modelo predijo correctamente que no era  clase B. Aquí es 0,9762, lo que significa que el modelo predijo correctamente 97,62%
#Pos Pred Value: este valor indica el porcentaje de veces que el modelo predijo  clase B correctamente cuando hizo esa predicción. En este caso 98,57% de las veces
#Neg Pred Value: mide cuántas veces el modelo predijo correctamente que no era  clase B cuando realmente no era  clase B. En este caso es 0,9535, lo que significa que el modelo predijo correctamente 95,35% de las veces que no era  clase B.



------------------------------------------------------------------------------------------




################################################
#Modelo RDA (Regularized Discriminant Analysis)#
################################################
library(klaR)

colnames(train_data)
train_data$Clase <- as.factor(train_data$Clase)

#Entrenar modelo RDA
rda_model <- klaR::rda(Clase ~ ., data = train_data)

#Predicción en el conjunto de prueba
rda_predictions <- predict(rda_model, test_data)$class

#Realizar predicciones sobre el conjunto de prueba
rda_predictions <- predict(rda_model, newdata = test_data)

#Obtener la predicción (predicciones de la clase)
predicted_classes <- rda_predictions$class
predicted_classes
length(predicted_classes)

#Obtener las verdaderas etiquetas (las clases reales en el conjunto de prueba)
true_classes <- as.factor(test_data$Clase)
true_classes
length(true_classes)


#Crear la matriz de confusión (trat predicho testing vs. trat real testing)
rda_conf_matrix <- confusionMatrix(predicted_classes, true_classes)
print(rda_conf_matrix)


#Comentar tabla segun salga: ejemplO:
#Accuracy: la exactitud global del modelo es del 95.58%, lo que significa que el modelo clasificó correctamente el 95.58% de las instancias en total. El intervalo de confianza del 95% para la exactitud está entre 89.98% y 98.55%, lo que indica un alto grado de confianza en la precisión de esta estimación.
#Sensitivity: la sensibilidad para la clase "B" (Benigno) es del 100%, lo que significa que el modelo identificó correctamente todas las instancias que eran realmente "B". No hubo falsos negativos para la clase "B".
#Specificity: la especificidad para la clase "B" es del 88.10%. Esto significa que el modelo predijo correctamente como "M" (Maligno) el 88.10% de las veces que la instancia no era "B". En otras palabras, de todas las instancias que eran realmente "M", el modelo las clasificó correctamente como "M" el 88.10% de las veces.
#Pos Pred Value: el valor predictivo positivo para la clase "B" es del 93.42%. Esto significa que, de todas las instancias que el modelo predijo como "B", el 93.42% realmente eran "B".
#Neg Pred Value: el valor predictivo negativo para la clase "B" es del 100%. Esto significa que, de todas las instancias que el modelo predijo como no "B" (es decir, como "M"), el 100% realmente eran "M".


------------------------------------------------------------------------------------------





##################################
#Modelo kNN (k-Nearest Neighbors)#
##################################

library(caret)


# Identificar variables con varianza cero
nzv <- nearZeroVar(train_data)

# Eliminar variables con varianza cero
train_data <- train_data[, -nzv]

#Crear un modelo de k-NN utilizando el paquete caret
knnModel <- train(Clase ~ .,
                  data = train_data,
                  method = "knn",
                  trControl = trainControl(method = "cv", number = 10),
                  preProcess = c("center", "scale"),
                  tuneLength = 30)
knnModel
#Corregir con el que salga: Accuracy was used to select the optimal model using the largest value. The final value used for the model was k = 9.


plot(knnModel)
#Comentar plot, ejemplo: A mayor número de vecinos peor va clasificando. En la gráfica vemos que k=9 es óptimo.

#Realizar predicciones en el conjunto de prueba utilizando el modelo entrenado
predictions <- predict(knnModel, newdata = test_data )
predictions

#Evaluar la precisión del modelo utilizando la matriz de confusión
knn_conf_matrix <- confusionMatrix(predictions, test_data$Clase)
print(knn_conf_matrix)

#Comentar tabla segun salga: ejemplO:
#Accuracy: la exactitud global del modelo es del 99.12%, lo que significa que el modelo clasificó correctamente el 99.12% de las instancias en total. El intervalo de confianza del 95% para la exactitud está entre 95.17% y 99.98%, lo que indica un alto grado de confianza en la precisión de esta estimación.
#Sensitivity: la sensibilidad para la clase "B" (Benigno) es del 100.00%. Esto significa que el modelo identificó correctamente el 100.00% de las instancias que eran realmente "B". No hubo ningún falso negativo (ninguna instancia que era "B" se predijo como "M").
#Specificity: la especificidad para la clase "B" es del 97.62%. Esto significa que el modelo predijo correctamente como "M" (Maligno) el 97.62% de las veces que la instancia no era "B". En otras palabras, de todas las instancias que eran realmente "M", el modelo las clasificó correctamente como "M" el 97.62% de las veces. Hubo falsos positivos (instancias que eran "M" pero se predijeron como "B").
#Pos Pred Value: el valor predictivo positivo para la clase "B" es del 98.61%. Esto significa que, de todas las instancias que el modelo predijo como "B", el 98.61% realmente eran "B".
#Neg Pred Value: el valor predictivo negativo para la clase "B" es del 100.00%. Esto significa que, de todas las instancias que el modelo predijo como no "B" (es decir, como "M"), el 100.00% realmente eran "M".

#Obtener probabilidades
probabilities_knn <- predict(knnModel, newdata = test_data, type = "prob")
probabilities_knn



------------------------------------------------------------------------------------------




##############################
#SVM (Support Vector Machine)#
#############################

library(caret)
  
  
#Modelo de SVM lineal utilizando el paquete caret
svmModelLineal <- train(Clase ~.,
                        data = train_data,
                        method = "svmLinear",
                        trControl = trainControl(method = "cv", number = 10),
                        preProcess = c("center", "scale"),
                        tuneGrid = expand.grid(C = seq(0, 2, length = 20)), 
                        prob.model = TRUE) 
svmModelLineal

  
plot(svmModelLineal) 
#Comentar plot ejemplo: Se ve como en el modelo de clasificacion, con valores mayores del cost va disminuyendo la precision, el primer valor es el que mejor funciona (C = 0.1052632).

  
#Realizar predicciones en el conjunto de prueba utilizando el modelo entrenado
predictions <- predict(svmModelLineal, newdata = test_data )
predictions
  
#Evaluar la precisión del modelo utilizando la matriz de confusión
smvModelLineal_conf_matrix <- confusionMatrix(predictions, test_data$Clase)
print(smvModelLineal_conf_matrix)

#Comentar tabla segun salga:   

#SVM lineal
probabilities_svm_linear <- predict(svmModelLineal, newdata = test_data, type = "prob")
probabilities_svm_linear



------------------------------------------------------------------------------------------




###################################
#Árbol de decisión (decision tree)#
###################################
  
library(caret)
library(rpart)
library(rattle)
  
# Crear un modelo de DT utilizando el paquete caret
 dtModel <- train(Clase ~.,
                  data = train_data,
                  method = "rpart",
                  trControl = trainControl(method = "cv", number = 10),
                  preProcess = c("center", "scale"),
                  tuneLength = 10)
dtModel
plot(dtModel)
  
fancyRpartPlot(dtModel$finalModel, type = 4, main = "Árbol de Decisión", sub="")
  
  
#Evaluar el modelo con el conjunto de prueba
predictions_raw <- predict(dtModel, newdata = test_data, type = "raw") # raw = clases
predictions_raw
  
  
#Evaluar la precisión del modelo utilizando la matriz de confusión
dtModel_conf_matrix <- confusionMatrix(predictions_raw, test_data$Clase)
print(dtModel_conf_matrix)


#Obtener probabilidades
probabilities_dt <- predict(dtModel, newdata = test_data, type = "prob")



------------------------------------------------------------------------------------------



############################
#Comparación de precisiones#
############################


#Función para extraer las métricas principales de la matriz de confusión
extract_metrics <- function(conf_matrix) {
  list(
    Accuracy = conf_matrix$overall["Accuracy"],
    Kappa = conf_matrix$overall["Kappa"],
    Sensitivity = conf_matrix$byClass["Sensitivity"],
    Specificity = conf_matrix$byClass["Specificity"],
    Pos_Pred_Value = conf_matrix$byClass["Pos Pred Value"],
    Neg_Pred_Value = conf_matrix$byClass["Neg Pred Value"]
  )
}

#Extraer las métricas de cada modelo
lda_metrics <- extract_metrics(lda_conf_matrix)
qda_metrics <- extract_metrics(qda_conf_matrix)
rda_metrics <- extract_metrics(rda_conf_matrix)
knn_metrics <- extract_metrics(knn_conf_matrix)
svm_lineal_metrics <- extract_metrics(smvModelLineal_conf_matrix)
dt_metrics <- extract_metrics(dtModel_conf_matrix)

#Dataframe comparativo
comparison_table <- data.frame(
  Model = c("LDA", "QDA", "RDA", "kNN", "SVM Lineal", "Árbol de decisión"),
  Accuracy = c(lda_metrics$Accuracy, qda_metrics$Accuracy, rda_metrics$Accuracy, 
               knn_metrics$Accuracy, svm_lineal_metrics$Accuracy, dt_metrics$Accuracy),
  Sensitivity = c(lda_metrics$Sensitivity, qda_metrics$Sensitivity, rda_metrics$Sensitivity,
                  knn_metrics$Sensitivity, svm_lineal_metrics$Sensitivity, dt_metrics$Sensitivity),
  Specificity = c(lda_metrics$Specificity, qda_metrics$Specificity, rda_metrics$Specificity,
                  knn_metrics$Specificity, svm_lineal_metrics$Specificity, dt_metrics$Specificity),
  Pos_Pred_Value = c(lda_metrics$Pos_Pred_Value, qda_metrics$Pos_Pred_Value, rda_metrics$Pos_Pred_Value,
                     knn_metrics$Pos_Pred_Value, svm_lineal_metrics$Pos_Pred_Value, dt_metrics$Pos_Pred_Value),
  Neg_Pred_Value = c(lda_metrics$Neg_Pred_Value, qda_metrics$Neg_Pred_Value, rda_metrics$Neg_Pred_Value,
                     knn_metrics$Neg_Pred_Value, svm_lineal_metrics$Neg_Pred_Value, dt_metrics$Neg_Pred_Value)
)

#Tabla comparativa
print(comparison_table)




