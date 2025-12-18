# Nextflow: Базовые кирпичики пайплайна

## Минимальный пайплайн - что нужно?

1. **Process** - аналог rule в Snakemake
2. **Channel** - поток данных между процессами
3. **Workflow** - собираем все вместе
4. **Params** - параметры запуска
5. **Config** - конфигурация окружения

## 1. Process - базовый кирпичик

```groovy
process FASTQC {
    // Входные данные
    input:
    tuple val(sample_id), path(reads)
    
    // Выходные данные
    output:
    path("*_fastqc.html")
    path("*_fastqc.zip")
    
    // Команда для выполнения
    script:
    """
    fastqc ${reads}
    """
}
```

### Сравнение с Snakemake

**Snakemake:**
```python
rule fastqc:
    input: "data/{sample}.fastq"
    output: "results/{sample}_fastqc.html"
    shell: "fastqc {input}"
```

**Nextflow:**
```groovy
process FASTQC {
    input: path(fastq)
    output: path("*_fastqc.html")
    script: "fastqc ${fastq}"
}
```

## 2. Input директивы

### Простой файл
```groovy
input:
path(reads)  // Один файл
```

### Tuple - несколько значений вместе
```groovy
input:
tuple val(sample_id), path(reads)
// sample_id = "SRR001", reads = "SRR001.fastq"
```

### Несколько входов
```groovy
input:
path(reads)
path(reference)
```

### Each - комбинации
```groovy
input:
each path(reference)  // Один reference для всех reads
path(reads)
```

## 3. Output директивы

### Простой вывод
```groovy
output:
path("result.txt")
```

### Tuple на выходе
```groovy
output:
tuple val(sample_id), path("${sample_id}.bam")
```

### Несколько выходов
```groovy
output:
path("*.html"), emit: html
path("*.zip"), emit: zip
```

`emit:` - даем имя выходу для использования дальше

## 4. Script директива

### Простой скрипт
```groovy
script:
"""
fastqc ${reads} -o results/
"""
```

### Многострочный
```groovy
script:
"""
mkdir -p results
fastqc ${reads} -o results/
multiqc results/ -o final/
"""
```

### С переменными
```groovy
script:
def prefix = sample_id
"""
spades.py -1 ${reads[0]} -2 ${reads[1]} -o ${prefix}_assembly
"""
```

## 5. Директивы процесса (Process directives)

### Conda окружение
```groovy
process PROKKA {
    conda 'bioconda::prokka=1.14.6'
    
    input:
    path(assembly)
    
    script:
    """
    prokka ${assembly}
    """
}
```

### Docker контейнер
```groovy
process SPADES {
    container 'staphb/spades:3.15.5'
    
    script:
    """
    spades.py -1 ${reads[0]} -2 ${reads[1]}
    """
}
```

### Ресурсы
```groovy
process ASSEMBLY {
    cpus 8
    memory '16 GB'
    time '2h'
    
    script:
    """
    spades.py -t ${task.cpus} -m ${task.memory.toGiga()}
    """
}
```

### Публикация результатов
```groovy
process RESULTS {
    publishDir "results/", mode: 'copy'
    
    input:
    path(files)
    
    output:
    path("*.html")
    
    script:
    """
    cp ${files} .
    """
}
```

## 6. Channels - создание

### Из файлов
```groovy
// Один файл
Channel.fromPath("data/sample.fastq")

// Все файлы по маске
Channel.fromPath("data/*.fastq")

// Парные риды
Channel.fromFilePairs("data/*_R{1,2}.fastq")
// Создаст: [sample_id, [R1.fastq, R2.fastq]]
```

### Из списка
```groovy
Channel.of('sample1', 'sample2', 'sample3')
```

### Из CSV файла
```groovy
Channel.fromPath('samples.csv')
    .splitCsv(header: true)
    .map { row -> tuple(row.sample_id, file(row.reads)) }
```

## 7. Workflow - собираем все вместе

### Базовая структура
```groovy
workflow {
    // 1. Создаем каналы
    reads_ch = Channel.fromFilePairs("data/*_R{1,2}.fastq")
    
    // 2. Запускаем процессы
    FASTQC(reads_ch)
    TRIMMING(reads_ch)
    ASSEMBLY(TRIMMING.out)
}
```

### Pipe оператор (|)
```groovy
workflow {
    Channel.fromFilePairs("data/*_R{1,2}.fastq") 
        | FASTQC 
        | TRIMMING 
        | ASSEMBLY
}
```

## 8. Params - параметры

```groovy
// В начале файла main.nf
params.reads = "data/*_R{1,2}.fastq"
params.outdir = "results"
params.threads = 4

workflow {
    reads_ch = Channel.fromFilePairs(params.reads)
    // ...
}
```

**Запуск с параметрами:**
```bash
nextflow run main.nf --reads "other_data/*.fastq" --threads 8
```

## 9. Минимальный полный пример

```groovy
#!/usr/bin/env nextflow

// Параметры
params.reads = "data/*_R{1,2}.fastq"
params.outdir = "results"

// Процесс 1: QC
process FASTQC {
    conda 'bioconda::fastqc=0.12.1'
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path("*_fastqc.{html,zip}")
    
    script:
    """
    fastqc ${reads}
    """
}

// Процесс 2: Trimming
process TRIMMOMATIC {
    conda 'bioconda::trimmomatic=0.39'
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R{1,2}.fastq")
    
    script:
    """
    trimmomatic PE ${reads[0]} ${reads[1]} \
        ${sample_id}_trimmed_R1.fastq ${sample_id}_unpaired_R1.fastq \
        ${sample_id}_trimmed_R2.fastq ${sample_id}_unpaired_R2.fastq \
        TRAILING:20 MINLEN:50
    """
}

// Workflow
workflow {
    // Создаем канал
    reads_ch = Channel.fromFilePairs(params.reads)
    
    // Запускаем процессы
    FASTQC(reads_ch)
    TRIMMOMATIC(reads_ch)
}
```

## 10. Config файл (nextflow.config)

```groovy
// nextflow.config
process {
    executor = 'local'
    cpus = 4
    memory = '8 GB'
    
    withName: ASSEMBLY {
        cpus = 16
        memory = '32 GB'
    }
}

conda {
    enabled = true
}

docker {
    enabled = false
}
```

## Запуск пайплайна

```bash
# Базовый запуск
nextflow run main.nf

# С параметрами
nextflow run main.nf --reads "data/*.fastq" --outdir "my_results"

# С профилем из конфига
nextflow run main.nf -profile docker

# Resume после ошибки
nextflow run main.nf -resume
```

## Структура проекта

```
my_pipeline/
├── main.nf              # Основной файл пайплайна
├── nextflow.config      # Конфигурация
├── modules/             # Модули (опционально)
│   ├── fastqc.nf
│   └── trimming.nf
├── bin/                 # Вспомогательные скрипты
│   └── helper.py
└── env/                 # Conda environment файлы
    └── fastqc.yml
```

## Чеклист минимального пайплайна

- ✅ Параметры через `params`
- ✅ Процессы с `input`, `output`, `script`
- ✅ Conda или Docker для зависимостей
- ✅ `publishDir` для результатов
- ✅ Workflow блок для связи процессов
- ✅ Config файл для настроек

## Готовы переписать Snakemake пайплайн? 🚀