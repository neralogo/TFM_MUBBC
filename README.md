# Trabajo de fin de máster 

## *Autor*: Nermina Logo Lendo

### *Máster Universitario en Bioinformática y Biología Computacional, UAM, 2022-2023*

Este repositorio contiene los códigos y manuales de las diferentes funciones utilizadas para el tabajo de fin de máster titulado: *Evaluación del rendimiento diagnóstico para la práctica clínica en base a datos genómicos.*

Este trabajo de fin de máster está comprendido dentro de un proyecto mayor que se lleva a cabo actualmente en el Hospital General Universitario Gregorio Marañón. Se trata, entonces, de una prueba piloto de este proyecto, por lo que no cuenta con toda la información que sería necesaria para poder realizar todos los análisis planteados inicialmente. 

En este trabajo se ha realizado un filtraje de variantes potencialmente patogénicas de novo de exoma y se ha aplicado un algoritmo de decisión tipo Random Forest para obtener variables de interés fenotípicas a la hora de tomar una decisión sobre la realización de la prueba genética de exoma. Los resultados indican que las variables de mayor interés para tomar la decisión sobre la realización del diagnóstico genético son: la edad de inicio de la marcha, edad de la madre, edad del padre y edad de inicio de los síntomas. Posteriormente se ha obtenido el árbol con menor error OOB y se ha realizado una clasificación con Bootstrap de los pacientes para obtener el porcentaje medio e intervalo de confianza de cada camino del árbol seleccionado.
