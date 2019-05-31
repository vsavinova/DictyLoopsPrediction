# DictyLoopsPrediction

## Данные
Геном взят из [dictyBase](http://dictybase.org/db/cgi-bin/dictyBase/download/download.pl?area=blast_databases&ID=dicty_chromosomal.gz)

## Tools
Использовался [bedtools](https://github.com/arq5x/bedtools2) makewindows с параметром -w 2000 для создания файла с геномными окнами размером в 2000 в формате .bed

## План проекта
1. С использованием литературных источников сформировать список кандидатов - свойств ДНК для предсказания
2. Исследовать и выбрать программы для подсчета свойств ДНК
3. Применить программы и подготовить список свойств 
4. Исследовать различные модели для задачи предсказания положения петель
5. Исследовать результаты моделей
6. Биологическая интерпретация результатов моделей

## Поиск мотивов
1. PWM из базы [CIS-BP](http://cisbp.ccbr.utoronto.ca/matchlist.php?versionNumber=1.02&orderby=1&familyTF=0&statusTF=0&studyID=0&dataSourceTF=0&speciesTF=Dictyostelium_discoideum), данные лежат в папке data/motifs/CisBP
2. PWM преобразована для autosome софта (часть "Преобразование PWM для autosome-софта" в notebook-е **src/motifs/motifs.ipynb**)
3. После преобразования нужно запустить скрипт **src/motifs/02_hit_motif.py** с параметрами:
    
    `--genome` (имя fasta файла с геномом без указания пути к файлу)
    
    `--genomeFolder` (папка, где лежит файл с геномом)
    
    `--background` (частоты встречаемости нуклеотидов в геноме)
    
    `--folderMotifs` (папка, где лежат PWMs для каждого мотива)
    
    `--folderThresholds` (папка для файлов с threshold-ми)
    
    `--pValue` (default 1e-5)
    
    Пример:
  ```bash
  python 02_hit_motif.py --genome dicty --genomeFolder ../../data/ --background 0.388,0.112,0.112,0.388 --folderMotifs ../../data/motifs/pwm/ --folderThresholds ../../data/motifs/motif_thresholds/ --outputFolder ../../data/motifs/tmp
  ```

4. В notebook-е **src/motifs/motifs.ipynb** считается обогащение позиций петель мотивами 

## G-квадруплексы
### G4Hunter
1. **G4-score**: в пределах от −4 до 4 (0 для A и T, >0 для G и <0 для C), вычисляется как сред арифметическое для score символов последовательности. G дает 1, если слева и справа от G A,T или C; GG - каждый G дает 2; >=4 подряд идущих G дает 4.
(Пример: score(atgaTTGGCGGGGAGAGGGAGGGGG) = (1+2*2+4*4+1+3*3+5*4 ) / 25 = 2.0)
2. **Запуск**: 
    ```bash 
    python G4Hunter.py -i <inputfile> -o <outputrepository> -w <window> -s <score threshold>
    ```
3. Подсчет G4, попавших в бины с петлями

## R-loops
### QmRLFS-finder
Запуск:
```bash 
    python QmRLFS-finder.py --log -bed -i ../data/dicty/dicty.fa -o ../data/R-loop/r_loops 
```