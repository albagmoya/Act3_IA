#las "-------------" separan los bloques de codigo para el rmd

rm(list=ls())

setwd("/Users/carme/OneDrive/Escritorio/R/Master_bioinformatica/Algoritmos/actividades/Act3_IA")

setwd("/home/albagmoya/Escriptori/IA/Act3_IA")

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

# Quitamos la columna Muestra
data_scaled <- data_scaled %>% dplyr::select(-Muestra)


#filtrado de datos columnas de varianza cero
filt_data_scaled <- data_scaled[,-nearZeroVar(data_scaled)]

------------------------------------------------------------------------------------------

# Dividir el conjunto de datos: entrenamiento (80%) y prueba (20%)

dataframe$Clase <- as.factor(dataframe$Clase)
levels(dataframe$Clase)  #Comprobar las clases existetes: "AGH" "CFB" "CGC" "CHC" "HPB"

set.seed(1995)

filt_data_scaled$Clase <- as.factor(filt_data_scaled$Clase)
train_index <- createDataPartition(filt_data_scaled$Clase, p = 0.8, list = FALSE)

train_data <- filt_data_scaled[train_index, ]
test_data <- filt_data_scaled[-train_index, ]

train_data$Clase <- factor(train_data$Clase, levels = c("AGH", "CFB", "CGC", "CHC", "HPB"))
test_data$Clase <- factor(test_data$Clase, levels = c("AGH", "CFB", "CGC", "CHC", "HPB"))

str(filt_data_scaled)


# Seleccionar las columnas numéricas para escalarlas
# numerical_columns_train <- train_data[, sapply(train_data, is.numeric)]
# scaled_train_data <- scale(numerical_columns_train)

# Convertir los datos escalados a un data.frame
# scaled_train_data <- as.data.frame(scaled_train_data)

# Unir los datos escalados con la variable Clase y Muestra en el conjunto de entrenamiento
# train_data <- cbind(scaled_train_data, Clase = train_data$Clase, Muestra = train_data$Muestra)

# Asegurarse de que el conjunto de datos de entrenamiento es un data.frame
# train_data <- as.data.frame(train_data)

# Escalar las columnas numéricas del conjunto de prueba
# numerical_columns_test <- test_data[, sapply(test_data, is.numeric)]
# scaled_test_data <- scale(numerical_columns_test)

# Convertir los datos escalados a un data.frame
# scaled_test_data <- as.data.frame(scaled_test_data)

# Unir los datos escalados con la variable Clase y Muestra en el conjunto de prueba
# test_data <- cbind(scaled_test_data, Clase = test_data$Clase, Muestra = test_data$Muestra)

# Asegurarse de que el conjunto de datos de prueba es un data.frame
# test_data <- as.data.frame(test_data)

# Confirmar la estructura de los datos
str(train_data)
str(test_data)

nearZeroVar(data_scaled, saveMetrics = TRUE)

colSums(apply(data_scaled, 2, var, na.rm = TRUE) == 0)


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
# La matriz de confusión muestra el número de predicciones correctas e incorrectas para cada clase:
#   
#   AGH: 29 casos fueron correctamente clasificados como AGH. No hubo errores.
# CFB: 60 casos fueron correctamente clasificados como CFB. Hubo 1 falso positivo (se predijo CFB cuando era CGC).
# CGC: 27 casos fueron correctamente clasificados como CGC. Hubo 1 falso negativo (se predijo HPB cuando era CGC).
# CHC: 27 casos fueron correctamente clasificados como CHC. No hubo errores.
# HPB: 14 casos fueron correctamente clasificados como HPB. No hubo errores.
# Estadísticas Globales
# Accuracy (Exactitud): La exactitud global del modelo es del 98.74%, lo que indica que clasificó correctamente el 98.74% de todos los casos. El intervalo de confianza del 95% (0.9553, 0.9985) sugiere una alta confianza en esta estimación.
# Kappa: El coeficiente Kappa de 0.9833 indica una excelente concordancia entre las predicciones del modelo y las etiquetas reales, más allá del azar.
# Valor P [Exactitud > Tasa de No Información]: El valor p < 2.2e-16 es extremadamente significativo. Esto indica que la exactitud del modelo es significativamente superior a la tasa de no información (0.3774).
# Estadísticas por Clase
# Sensibilidad (Recall):
#   AGH: 100% (todos los casos AGH reales fueron identificados correctamente).
# CFB: 100% (todos los casos CFB reales fueron identificados correctamente).
# CGC: 96.43% (casi todos los casos CGC reales fueron identificados correctamente).
# CHC: 100% (todos los casos CHC reales fueron identificados correctamente).
# HPB: 93.33% (la mayoría de los casos HPB reales fueron identificados correctamente).
# Especificidad:
#   AGH: 100% (todos los casos que no eran AGH fueron clasificados correctamente como no AGH).
# CFB: 98.99% (casi todos los casos que no eran CFB fueron clasificados correctamente como no CFB).
# CGC: 99.24% (casi todos los casos que no eran CGC fueron clasificados correctamente como no CGC).
# CHC: 100% (todos los casos que no eran CHC fueron clasificados correctamente como no CHC).
# HPB: 100% (todos los casos que no eran HPB fueron clasificados correctamente como no HPB).
# Valor Predictivo Positivo (PPV):
#   AGH: 100% (todos los casos predichos como AGH realmente eran AGH).
# CFB: 98.36% (casi todos los casos predichos como CFB realmente eran CFB).
# CGC: 96.43% (casi todos los casos predichos como CGC realmente eran CGC).
# CHC: 100% (todos los casos predichos como CHC realmente eran CHC).
# HPB: 100% (todos los casos predichos como HPB realmente eran HPB).
# Valor Predictivo Negativo (NPV):
#   AGH: 100% (todos los casos predichos como no AGH realmente no eran AGH).
# CFB: 100% (todos los casos predichos como no CFB realmente no eran CFB).
# CGC: 99.24% (casi todos los casos predichos como no CGC realmente no eran CGC).
# CHC: 100% (todos los casos predichos como no CHC realmente no eran CHC).
# HPB: 99.31% (casi todos los casos predichos como no HPB realmente no eran HPB).


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

train_data <- train_data[, sapply(train_data, is.numeric) | names(train_data) == "Clase"]

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

set.seed(1995)
#Crear un modelo de k-NN utilizando el paquete caret
knnModel <- train(Clase ~ .,
                  data = train_data,
                  method = "knn",
                  trControl = trainControl(method = "cv", number = 10),
                  preProcess = c("center", "scale"),
                  tuneLength = 30)
knnModel
# Accuracy was used to select the optimal model using the largest value. The final value used for the model was k = 15.


plot(knnModel)
#A mayor número de vecinos peor va clasificando. En la gráfica vemos que k=15 es óptimo.

#Realizar predicciones en el conjunto de prueba utilizando el modelo entrenado
predictions <- predict(knnModel, newdata = test_data )
predictions

#Evaluar la precisión del modelo utilizando la matriz de confusión
knn_conf_matrix <- confusionMatrix(predictions, test_data$Clase)
print(knn_conf_matrix)

#Tabla comentada:

# La matriz de confusión muestra el número de predicciones correctas e incorrectas para cada clase:

# AGH: 28 casos fueron correctamente clasificadas como AGH. Hubo 1 falso negativo (se predijo CHC cuando era AGH).
# CFB: 60 casos fueron correctamente clasificadas como CFB. No hubo errores.
# CGC: 28 casos fueron correctamente clasificadas como CGC. Hubo 1 falso positivo (se predijo CGC cuando era HPB).
# CHC: 27 casos fueron correctamente clasificadas como CHC. Hubo 1 falso negativo (se predijo AGH cuando era CHC).
# HPB: 14 casos fueron correctamente clasificadas como HPB. No hubo errores.
# Estadísticas Globales
# Accuracy (Exactitud): La exactitud global del modelo es del 98.74%, lo que indica que clasificó correctamente el 98.74% de todos los casos. El intervalo de confianza del 95% (0.9553, 0.9985) sugiere una alta confianza en esta estimación.
# Kappa: El coeficiente Kappa de 0.9833 indica una excelente concordancia entre las predicciones del modelo y las etiquetas reales, más allá del azar.
# Valor P [Exactitud > Tasa de No Información]: El valor p < 2.2e-16 es extremadamente significativo. Esto indica que la exactitud del modelo es significativamente superior a la tasa de no información (0.3774), que representa la exactitud que se obtendría si siempre se predijera la clase más frecuente.
# Estadísticas por Clase
# Sensibilidad (Recall):
#   AGH: 96.55% (casi todos los casos AGH reales fueron identificadas correctamente).
# CFB: 100% (todos los casos CFB reales fueron identificadas correctamente).
# CGC: 100% (todos los casos CGC reales fueron identificadas correctamente).
# CHC: 100% (todos los casos CHC reales fueron identificadas correctamente).
# HPB: 93.33% (la mayoría de los casos HPB reales fueron identificadas correctamente).
# Especificidad:
#   AGH: 100% (todos los casos que no eran AGH fueron clasificadas correctamente como no AGH).
# CFB: 100% (todos los casos que no eran CFB fueron clasificadas correctamente como no CFB).
# CGC: 99.24% (casi todos los casos que no eran CGC fueron clasificadas correctamente como no CGC).
# CHC: 99.24% (casi todos los casos que no eran CHC fueron clasificadas correctamente como no CHC).
# HPB: 100% (todos los casos que no eran HPB fueron clasificadas correctamente como no HPB).
# Valor Predictivo Positivo (PPV):
#   AGH: 100% (todos los casos predichas como AGH realmente eran AGH).
# CFB: 100% (todos los casos predichas como CFB realmente eran CFB).
# CGC: 96.55% (la mayoría de las casos predichas como CGC realmente eran CGC).
# CHC: 96.43% (la mayoría de las casos predichas como CHC realmente eran CHC).
# HPB: 100% (todos los casos predichas como HPB realmente eran HPB).
# Valor Predictivo Negativo (NPV):
#   AGH: 99.24% (casi todos los casos predichas como no AGH realmente no eran AGH).
# CFB: 100% (todos los casos predichas como no CFB realmente no eran CFB).
# CGC: 100% (todos los casos predichas como no CGC realmente no eran CGC).
# CHC: 100% (todos los casos predichas como no CHC realmente no eran CHC).
# HPB: 99.31% (casi todos los casos predichas como no HPB realmente no eran HPB).



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

#En el gráfico se observa que el valor del cost mejor está cercano a 1.5 (C = 1.368421).

  
#Realizar predicciones en el conjunto de prueba utilizando el modelo entrenado
predictions <- predict(svmModelLineal, newdata = test_data )
predictions
  
#Evaluar la precisión del modelo utilizando la matriz de confusión
smvModelLineal_conf_matrix <- confusionMatrix(predictions, test_data$Clase)
print(smvModelLineal_conf_matrix)

# La matriz de confusión muestra el número de predicciones correctas e incorrectas para cada clase:

# AGH: 29 casos fueron correctamente clasificados como AGH. No hubo errores.
# CFB: 60 casos fueron correctamente clasificados como CFB. Hubo 2 falsos positivos (se predijo CFB cuando era CGC).
# CGC: 28 casos fueron correctamente clasificados como CGC. Hubo 16 falsos positivos (se predijo CGC cuando era CHC) y 8 falsos positivos (se predijo CGC cuando era HPB).
# CHC: 9 casos fueron correctamente clasificados como CHC.
# HPB: 7 casos fueron correctamente clasificados como HPB.
# Estadísticas Globales
# Accuracy (Exactitud): La exactitud global del modelo es del 83.65%, lo que indica que clasificó correctamente el 83.65% de todos los casos. El intervalo de confianza del 95% (0.7697, 0.8903) sugiere una confianza moderada en esta estimación. Notablemente menor que el modelo anterior.
# Kappa: El coeficiente Kappa de 0.7815 indica una buena concordancia entre las predicciones del modelo y las etiquetas reales, aunque no tan alta como el modelo anterior.
# Valor P [Exactitud > Tasa de No Información]: El valor p < 2.2e-16 es extremadamente significativo. Esto indica que la exactitud del modelo es significativamente superior a la tasa de no información (0.3774).
# Estadísticas por Clase
# Sensibilidad (Recall):
# AGH: 100% (todos los casos AGH reales fueron identificados correctamente).
# CFB: 100% (todos los casos CFB reales fueron identificados correctamente).
# CGC: 100% (todos los casos CGC reales fueron identificados correctamente).
# CHC: 33.33% (solo un tercio de los casos CHC reales fueron identificados correctamente).
# HPB: 46.67% (menos de la mitad de los casos HPB reales fueron identificados correctamente).
# Especificidad:
# AGH: 100% (todos los casos que no eran AGH fueron clasificados correctamente como no AGH).
# CFB: 97.98% (casi todos los casos que no eran CFB fueron clasificados correctamente como no CFB).
# CGC: 81.68% (un número considerable de casos que no eran CGC fueron incorrectamente clasificados como CGC).
# CHC: 100% (todos los casos que no eran CHC fueron clasificados correctamente como no CHC).
# HPB: 100% (todos los casos que no eran HPB fueron clasificados correctamente como no HPB).
# Valor Predictivo Positivo (PPV):
# AGH: 100% (todos los casos predichos como AGH realmente eran AGH).
# CFB: 96.77% (casi todos los casos predichos como CFB realmente eran CFB).
# CGC: 53.85% (algo más de la mitad de los casos predichos como CGC realmente eran CGC).
# CHC: 100% (todos los casos predichos como CHC realmente eran CHC).
# HPB: 100% (todos los casos predichos como HPB realmente eran HPB).
# Valor Predictivo Negativo (NPV):
# AGH: 100% (todos los casos predichos como no AGH realmente no eran AGH).
# CFB: 100% (todos los casos predichos como no CFB realmente no eran CFB).
# CGC: 100% (todos los casos predichos como no CGC realmente no eran CGC).
# CHC: 88% (la mayoría de los casos predichos como no CHC realmente no eran CHC).
# HPB: 94.74% (casi todos los casos predichos como no HPB realmente no eran HPB).
# Prevalencia, Tasa de Detección, Prevalencia de Detección, Exactitud Balanceada: (Se interpretan de la misma manera que en el ejemplo anterior, pero con los valores correspondientes).

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

str(train_data)

any(is.na(train_data))

formula <- as.formula("Clase ~ .")
all.vars(formula) %in% colnames(train_data)

missing_vars <- setdiff(all.vars(formula), colnames(train_data))
print(missing_vars)

predictors <- setdiff(colnames(train_data), "Clase")
formula <- as.formula(paste("Clase ~", paste(predictors, collapse = " + ")))





predictors <- setdiff(colnames(train_data), "Clase")
print(predictors)

formula <- as.formula(paste("Clase ~", paste(predictors, collapse = " + ")))
print(formula)

dtModel <- train(
  formula, 
  data = train_data,
  method = "rpart",
  trControl = trainControl(method = "cv", number = 10),
  preProcess = c("center", "scale"),
  tuneLength = 10
)


  
  
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




