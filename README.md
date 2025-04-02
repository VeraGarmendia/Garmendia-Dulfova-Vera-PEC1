# README
Esta PEC tiene como objetivo planificar y ejecutar una versión simplificada del proceso de análisis de datos
ómicos en metabolómica, utilizando R y sus librerias, como Bioconductor, y creando un objeto SummarizedExperiment.

Los datos consisten en tres archivos separados:
1. Feature data: 1541 variables o metabolitos, 45 muestras.
2. Metadata: 45 filas, dos columnas (nombre de muestra, nombre de tratamiento).
3. Metabolites name: 1541 filas, 3 columnas (nombre original del metabolito, ID de PubChem y ID de
KEGG).

Estos datos se han obtenido del repositorio Metabolomics Workbench con ID ST000291 (“LC-MS Based Approaches to Investigate Metabolomic Differences in the
Urine of Young Women after Drinking Cranberry Juice or Apple Juice”).
