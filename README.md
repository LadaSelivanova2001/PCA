# PCA
Создан класс "PCA".
В нем, использована математическая библиотека, разработанная мной, а также реализованы методы:

  * Предобработки данных (центрирование и шкалирование). Для шкалирования использовать формулу, где деление производится на n-1.
  * Разложения матрицы алгоритмом NIPALS разложения матрицы. Вход - матрица X, выход - матрица счетов T (scores) , матрица весов P (loadings), матрица остатков E
  * Вычисления размахов и отклонений - считаются для каждого образца (строки в матрице).
  * Вычисления полной и объясненной дисперсий остатков в обучении.
  
  Файлы loadings.txt и scores.txt созданы для проверки расчетов.
