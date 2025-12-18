# Nextflow: Dataflow vs Graph-based

## Две философии пайплайнов

### Snakemake - Graph-based (Pull model)
**Логика:** "Мне нужен финальный файл, какие файлы нужны для его создания?"

```
Финальный файл ← Промежуточный файл ← Входной файл
     ↑                    ↑                   ↑
   Правило 2          Правило 1         Существует
```

Snakemake идет **от результата к началу** (backward chaining)

### Nextflow - Dataflow (Push model)
**Логика:** "Есть данные, что с ними делать дальше?"

```
Входные данные → Процесс 1 → Процесс 2 → Результат
                     ↓            ↓           ↓
                  Channel     Channel     Channel
```

Nextflow идет **от начала к результату** (forward chaining)

## Концепция динамического программирования

### Что это значит в контексте Nextflow?

**Динамическое программирование** = решения принимаются на ходу, основываясь на текущих данных

```groovy
// Nextflow - данные текут по каналам
process FASTQC {
    input:
    tuple val(sample_id), path(reads)  // Что придет - то обработаем
    
    output:
    path("*_fastqc.html")
    
    script:
    """
    fastqc ${reads}
    """
}

// Процесс НЕ ЗНАЕТ заранее сколько образцов придет!
// Он запустится столько раз, сколько элементов в канале
```

### Snakemake - статическая природа

```python
# Snakemake - нужно заранее знать все файлы
rule all:
    input:
        expand("results/{sample}_fastqc.html", 
               sample=["S1", "S2", "S3"])  # Жестко задано!
```

## Ключевые отличия

| Аспект | Snakemake | Nextflow |
|--------|-----------|----------|
| Парадигма | Pull (от результата) | Push (от данных) |
| Планирование | До запуска (DAG строится заранее) | Во время выполнения |
| Входные файлы | Должны быть известны | Могут определяться динамически |
| Параллелизм | По wildcards | По элементам в каналах |
| Язык | Python-based | Groovy-based |

## Визуализация: как это работает

### Snakemake
```
1. Определяем финальные файлы в rule all
2. Строим полный граф зависимостей (DAG)
3. Запускаем задачи согласно графу
```

```python
rule all:
    input: "final.txt"  # Хочу этот файл!

# Snakemake строит граф:
# final.txt ← process.txt ← raw.txt
# Запускает: правило для raw.txt → правило для process.txt → правило для final.txt
```

### Nextflow
```
1. Создаем каналы с данными
2. Данные текут через процессы
3. Каждый процесс запускается, когда получает данные
```

```groovy
// Создаем канал с данными
Channel.fromPath("data/*.fastq") | FASTQC | MULTIQC

// Данные текут:
// files → FASTQC → results → MULTIQC → report
// Каждый файл обрабатывается независимо и параллельно
```

## Пример: условная логика

### Snakemake (сложнее)
```python
def get_input(wildcards):
    if config["has_assembly"]:
        return f"assembly/{wildcards.sample}.fasta"
    else:
        return f"reads/{wildcards.sample}.fastq"

rule process:
    input: get_input
    output: "results/{sample}.txt"
    shell: "process.sh {input} > {output}"
```

### Nextflow (естественнее)
```groovy
reads_ch = Channel.fromPath("reads/*.fastq")
assembly_ch = Channel.fromPath("assembly/*.fasta")

// Просто объединяем каналы
input_ch = reads_ch.mix(assembly_ch)

process PROCESS {
    input:
    path(input_file)
    
    output:
    path("result.txt")
    
    script:
    """
    process.sh ${input_file} > result.txt
    """
}

input_ch | PROCESS
```

## Когда использовать что?

### Snakemake лучше если:
- Четко определенная структура проекта
- Файлы известны заранее
- Нужна интеграция с Python экосистемой
- Простой проект, понятная структура

### Nextflow лучше если:
- Динамическое количество образцов
- Сложная условная логика
- Масштабирование на кластеры/облака
- Нужна гибкость потоков данных

## Практическая аналогия

**Snakemake** = Рецепт блюда
- Знаешь ингредиенты → знаешь все шаги → готовишь

**Nextflow** = Конвейер на производстве  
- Детали идут по ленте → каждая станция делает свою работу → готово

## Главное

- **Snakemake** - думает о **файлах** (graph of files)
- **Nextflow** - думает о **данных** (flow of data)

Обе системы мощные, выбор зависит от задачи!