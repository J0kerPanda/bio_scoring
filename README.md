# Утилита командной строки для сравнения нуклеотидных/белковых последовательностей на языке Scala

1. Инструкция по сборке:  
  1.1. Установить Scala Build Tool (sbt):
    ```
    echo "deb https://dl.bintray.com/sbt/debian /" | sudo tee -a /etc/apt/sources.list.d/sbt.list
    sudo apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 2EE0EA64E40A89B84B2DF73499E82A75642AC823
    sudo apt-get update
    sudo apt-get install sbt
    ```  
    1.2. Запустить команду ```sbt assembly``` в корне проекта. Результатом её работы является файл с расширением ```.jar``` в директории ```target/scala-2.12```


2. Инструкция по запуску:  
  Для запуска проекта используется команда ```java -jar <file> [options]```. Предполагается, что на компьютере установлен JDK версии не ниже 8.
  Дополнительные опции запуска:
    ```
    --seq1 <value>           Путь к файлу первой последовательности
    --seq2 <value>           Путь к файлу второй последовательности
    -s, --sequenceType <value>
                             Тип последовательностей (n - нуклеотиды, p - протеины)
    -g, --gap <value>        Штраф за пропуски
    -e, --extended_gap <value> 
                             Штраф за последовательные пропуски(defaults to gap)
    --group <value>          Размер групп, на которые разделяются выровненные последовательности
    --help                   Напечатать вспомогательное сообщение
    ```
  
3. Пример работы алгоритма: 

Строка 1: ```ATTAAAT```

Строка 2: ```AAAT```  
Тип последовательностей: белки  
Матрица: DNAFull  
Штраф за пропуски: -5 
Штраф за последовательные пропуски: -5

Целью алгоритма является нахождение выравнивания, при котором результат (```score```), соответствующий разнице между строками, будет минимальным. Данный параметр рассчитывается исходя из того, что любое промежуточное выравнивание можно получить, добавляя последующую букву или пропуск в каждую строчку другого выравнивания за исключением случая добавления пропусков в каждую из строчек. Базовыми промежуточными выравниваниями считаются выравнивания, состоящие только из пропусков.  

Поскольку получение нового выравнивания из предыдущих неоднозначно, выбор происходит на основе максимизации параметра ```score```. Добавление пропуска в том или ином направлении в предыдущее выравнивание добавляет к параметру штраф за пропуск, а изменение параметра при фиксации двух символов происходит на основе весовой матрицы.
  
  Базовый шаг - заполнение базовых случаев (вставка пропусков): 
  
  |   	|     	| A  	| A   	| A   	| T   	|
  |---	|-----	|----	|-----	|-----	|-----	|
  |   	| 0   	| -5 	| -10 	| -15 	| -20 	|
  | A 	| -5  	|    	|     	|     	|     	|
  | T 	| -10 	|    	|     	|     	|     	|
  | T 	| -15 	|    	|     	|     	|     	|
  | A 	| -20 	|    	|     	|     	|     	|
  | A 	| -25 	|    	|     	|     	|     	|
  | A 	| -30 	|    	|     	|     	|     	|
  | T 	| -35 	|    	|     	|     	|     	|
  
  Шаг 1 - заполнение промежуточных выравниваний для 1-ых букв:  
  
  |   	|     	| A   	| A   	| A   	| T   	|
  |---	|-----	|-----	|-----	|-----	|-----	|
  |   	| 0   	| -5  	| -10 	| -15 	| -20 	|
  | A 	| -5  	|  5  	| 0   	| -5  	| -10 	|
  | T 	| -10 	| 0   	|     	|     	|     	|
  | T 	| -15 	| -5  	|     	|     	|     	|
  | A 	| -20 	| -10 	|     	|     	|     	|
  | A 	| -25 	| -15 	|     	|     	|     	|
  | A 	| -30 	| -20 	|     	|     	|     	|
  | T 	| -35 	| -25 	|     	|     	|     	|
  
  Шаг 2 - заполнение промежуточных выравниваний для 2-ых букв:  

  |   	|     	| A   	| A   	| A   	| T   	|
  |---	|-----	|-----	|-----	|-----	|-----	|
  |   	| 0   	| -5  	| -10 	| -15 	| -20 	|
  | A 	| -5  	|  5  	| 0   	| -5  	| -10 	|
  | T 	| -10 	| 0   	| 1   	| -4  	| 0   	|
  | T 	| -15 	| -5  	| -4  	|     	|     	|
  | A 	| -20 	| -10 	| 0   	|     	|     	|
  | A 	| -25 	| -15 	| -5  	|     	|     	|
  | A 	| -30 	| -20 	| -10 	|     	|     	|
  | T 	| -35 	| -25 	| -15 	|     	|     	|
  
  Шаг 3 - заполнение промежуточных выравниваний для 3-их букв:  

  |   	|     	| A   	| A   	| A   	| T   	|
  |---	|-----	|-----	|-----	|-----	|-----	|
  |   	| 0   	| -5  	| -10 	| -15 	| -20 	|
  | A 	| -5  	|  5  	| 0   	| -5  	| -10 	|
  | T 	| -10 	| 0   	| 1   	| -4  	| 0   	|
  | T 	| -15 	| -5  	| -4  	| -3  	| 1   	|
  | A 	| -20 	| -10 	| 0   	| 1   	|     	|
  | A 	| -25 	| -15 	| -5  	| 5   	|     	|
  | A 	| -30 	| -20 	| -10 	| 0   	|     	|
  | T 	| -35 	| -25 	| -15 	| -5  	|     	|
  
  Шаг 4 - заполнение промежуточных выравниваний для 4-ых букв: 
  
  |   	|     	| A   	| A   	| A   	| T   	|
  |---	|-----	|-----	|-----	|-----	|-----	|
  |   	| 0   	| -5  	| -10 	| -15 	| -20 	|
  | A 	| -5  	| 5   	| 0   	| -5  	| -10 	|
  | T 	| -10 	| 0   	| 1   	| -4  	| 0   	|
  | T 	| -15 	| -5  	| -4  	| -3  	| 1   	|
  | A 	| -20 	| -10 	| 0   	| 1   	| -4  	|
  | A 	| -25 	| -15 	| -5  	| 5   	| 0   	|
  | A 	| -30 	| -20 	| -10 	| 0   	| 1   	|
  | T 	| -35 	| -25 	| -15 	| -5  	| 5   	|
  
  Таким образом, в результате работы алгоритма параметр ```score``` определён как 5. Для выведения соответствующего выравнивания достаточно запоминать индексы колонок и строк, из которых осуществляется переход, при получении каждого нового выравнивания:
  s1: ATTAAAT
  s2: A--AA-T
 
