# Nextflow Channels: Потоки данных

## Что такое Channel?

**Channel** = труба, по которой текут данные между процессами

```
Файлы → Channel → Process 1 → Channel → Process 2 → Результат
```

Каналы могут быть:
- **Value channels** - один элемент, используется много раз
- **Queue channels** - поток элементов, каждый используется один раз

## 1. Создание каналов

### fromPath - файлы по маске

```groovy
// Один файл
Channel.fromPath('data/sample.fastq')

// Множество файлов
Channel.fromPath('data/*.fastq')

// Рекурсивно
Channel.fromPath('data/**/*.fastq')

// С сортировкой
Channel.fromPath('data/*.fastq').toSortedList()
```

### fromFilePairs - парные файлы

```groovy
// Автоматически находит пары
Channel.fromFilePairs('data/*_R{1,2}.fastq')
// Результат: [sample_id, [R1_file, R2_file]]

// Пример:
// data/SRR001_R1.fastq, data/SRR001_R2.fastq
// → ['SRR001', [SRR001_R1.fastq, SRR001_R2.fastq]]
```

### of - из значений

```groovy
Channel.of(1, 2, 3, 4, 5)

Channel.of(['sample1', 'path/to/file1.fq'],
           ['sample2', 'path/to/file2.fq'])
```

### fromList - из списка

```groovy
def samples = ['sample1', 'sample2', 'sample3']
Channel.fromList(samples)
```

### from - универсальный (deprecated, но полезно знать)

```groovy
Channel.from('sample1', 'sample2', 'sample3')
```

## 2. Операторы каналов

### map - трансформация

```groovy
Channel.fromPath('data/*.fastq')
    .map { file -> [file.simpleName, file] }
// file.simpleName - имя файла без расширения
// Результат: ['SRR001', SRR001.fastq]
```

### filter - фильтрация

```groovy
Channel.fromPath('data/*.fastq')
    .filter { it.size() > 1000000 }  // Только файлы > 1MB
```

### unique - уникальные значения

```groovy
Channel.of('A', 'B', 'A', 'C', 'B')
    .unique()
// Результат: A, B, C
```

### flatten - развернуть вложенные списки

```groovy
Channel.of([1, 2], [3, 4], [5])
    .flatten()
// Результат: 1, 2, 3, 4, 5
```

### collect - собрать все в список

```groovy
Channel.of(1, 2, 3, 4, 5)
    .collect()
// Результат: [1, 2, 3, 4, 5]
```

### splitCsv - парсинг CSV

```groovy
Channel.fromPath('samples.csv')
    .splitCsv(header: true, sep: ',')
    .map { row -> tuple(row.sample_id, file(row.reads)) }
```

**samples.csv:**
```
sample_id,reads
SRR001,data/SRR001.fastq
SRR002,data/SRR002.fastq
```

## 3. Комбинирование каналов

### mix - объединить каналы

```groovy
reads_ch = Channel.fromPath('reads/*.fastq')
assembly_ch = Channel.fromPath('assembly/*.fasta')

combined_ch = reads_ch.mix(assembly_ch)
// Все файлы в одном канале
```

### join - объединить по ключу

```groovy
reads_ch = Channel.of(['sample1', 'reads.fq'],
                       ['sample2', 'reads.fq'])

qc_ch = Channel.of(['sample1', 'qc.html'],
                    ['sample2', 'qc.html'])

reads_ch.join(qc_ch)
// Результат: ['sample1', 'reads.fq', 'qc.html']
//           ['sample2', 'reads.fq', 'qc.html']
```

### combine - декартово произведение

```groovy
samples = Channel.of('sample1', 'sample2')
refs = Channel.of('ref1.fa', 'ref2.fa')

samples.combine(refs)
// Результат: ['sample1', 'ref1.fa']
//           ['sample1', 'ref2.fa']
//           ['sample2', 'ref1.fa']
//           ['sample2', 'ref2.fa']
```

## 4. Многократное использование

### Value channel (one to many)

```groovy
// Один reference для всех
reference_ch = Channel.value(file('reference.fa'))

workflow {
    reads_ch = Channel.fromPath('*.fastq')
    
    // reference будет использоваться для каждого reads
    ALIGN(reads_ch, reference_ch)
}
```

## 5. Практические задачки 🎯

### Задачка 1: Простая трансформация

**Дано:** Файлы `sample1.fastq`, `sample2.fastq`, `sample3.fastq`

**Задача:** Создайте канал, который для каждого файла создаст tuple `[sample_name, file_path]`

```groovy
// ??? ВАШ КОД
Channel.fromPath('*.fastq')
    .map { ??? }
    .view()  // Показать результат

// Ожидаемый результат:
// ['sample1', sample1.fastq]
// ['sample2', sample2.fastq]
// ['sample3', sample3.fastq]
```

<details>
<summary>Решение</summary>

```groovy
Channel.fromPath('*.fastq')
    .map { file -> [file.simpleName, file] }
    .view()
```
</details>

---

### Задачка 2: Парсинг CSV

**Дано:** `samples.csv`
```csv
sample,forward,reverse,group
SRR001,SRR001_R1.fq,SRR001_R2.fq,control
SRR002,SRR002_R1.fq,SRR002_R2.fq,treatment
```

**Задача:** Создайте канал формата `[sample, group, [forward, reverse]]`

```groovy
// ??? ВАШ КОД
Channel.fromPath('samples.csv')
    .splitCsv(???)
    .map { ??? }
    .view()
```

<details>
<summary>Решение</summary>

```groovy
Channel.fromPath('samples.csv')
    .splitCsv(header: true)
    .map { row -> 
        [row.sample, 
         row.group, 
         [file(row.forward), file(row.reverse)]]
    }
    .view()
```
</details>

---

### Задачка 3: Фильтрация

**Дано:** Канал с файлами разного размера

**Задача:** Оставьте только файлы больше 1MB и покажите их размер в GB

```groovy
Channel.fromPath('data/*.fastq')
    .??? // фильтр по размеру > 1MB
    .map { ??? } // [filename, size_in_GB]
    .view()
```

<details>
<summary>Решение</summary>

```groovy
Channel.fromPath('data/*.fastq')
    .filter { it.size() > 1_000_000 }
    .map { file -> 
        [file.name, file.size() / 1_000_000_000]
    }
    .view()
```
</details>

---

### Задачка 4: Комбинирование каналов

**Дано:** 
- Канал с образцами: `['sample1'], ['sample2']`
- Канал с референсами: `['hg38.fa'], ['mm10.fa']`

**Задача:** Создайте все комбинации `[sample, reference]`

```groovy
samples_ch = Channel.of(['sample1'], ['sample2'])
refs_ch = Channel.of(['hg38.fa'], ['mm10.fa'])

// ??? ВАШ КОД
samples_ch.???
```

<details>
<summary>Решение</summary>

```groovy
samples_ch = Channel.of(['sample1'], ['sample2'])
refs_ch = Channel.of(['hg38.fa'], ['mm10.fa'])

samples_ch
    .combine(refs_ch)
    .view()

// Результат:
// ['sample1', 'hg38.fa']
// ['sample1', 'mm10.fa']
// ['sample2', 'hg38.fa']
// ['sample2', 'mm10.fa']
```
</details>

---

### Задачка 5: Группировка по ключу

**Дано:** Канал с результатами от разных образцов в разных репликах
```groovy
results = Channel.of(
    ['sample1', 'rep1', 'file1.txt'],
    ['sample1', 'rep2', 'file2.txt'],
    ['sample2', 'rep1', 'file3.txt'],
    ['sample2', 'rep2', 'file4.txt']
)
```

**Задача:** Сгруппируйте файлы по sample_id

```groovy
results
    .map { sample, rep, file -> [sample, file] }
    .groupTuple()
    .view()

// Ожидается:
// ['sample1', [file1.txt, file2.txt]]
// ['sample2', [file3.txt, file4.txt]]
```

## 6. Отладка каналов

### view() - показать содержимое

```groovy
Channel.fromPath('*.fastq')
    .view()  // Печатает каждый элемент
```

### view() с трансформацией

```groovy
Channel.fromPath('*.fastq')
    .view { "Processing file: $it" }
```

### ifEmpty() - проверка на пустоту

```groovy
Channel.fromPath('data/*.fastq')
    .ifEmpty { error "No FASTQ files found!" }
```

### count() - подсчет элементов

```groovy
Channel.fromPath('*.fastq')
    .count()
    .view { "Total files: $it" }
```

## 7. Частые паттерны

### Чтение CSV и создание tuple

```groovy
def samples_ch = Channel
    .fromPath(params.input_csv)
    .splitCsv(header: true)
    .map { row -> 
        tuple(
            row.sample_id,
            file(row.read1),
            file(row.read2)
        )
    }
```

### One reference, many samples

```groovy
reference = Channel.value(file(params.reference))
samples = Channel.fromPath('*.fastq')

process ALIGN {
    input:
    path(reads)
    path(ref)  // Будет одинаковый для всех
    
    script:
    "bwa mem $ref $reads"
}

workflow {
    ALIGN(samples, reference)
}
```

### Условное разделение потока

```groovy
Channel.fromPath('*.fastq')
    .branch {
        large: it.size() > 1_000_000
        small: it.size() <= 1_000_000
    }
    .set { result }

result.large | PROCESS_LARGE
result.small | PROCESS_SMALL
```

## Главное про каналы

- ✅ Каналы = потоки данных между процессами
- ✅ Queue channels используются один раз
- ✅ Value channels можно переиспользовать
- ✅ Операторы map, filter, join для трансформации
- ✅ view() для отладки
- ✅ CSV файлы - стандартный способ задать образцы

## Следующий шаг

Попробуйте переписать ваш Snakemake пайплайн, используя:
1. Process для каждого правила
2. Channel.fromPath или splitCsv для входных данных  
3. Операторы каналов для трансформации
4. Workflow для связи процессов

Удачи! 🚀