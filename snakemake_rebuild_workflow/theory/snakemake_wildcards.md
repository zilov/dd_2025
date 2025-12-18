# Snakemake Wildcards

## Зачем нужны wildcards?

**Проблема:** У нас 10 образцов, и мы не хотим писать 10 одинаковых правил

**Решение:** Wildcards - это переменные в путях файлов, которые Snakemake автоматически подставляет

## Базовый синтаксис

```python
rule trim_reads:
    input:
        r1 = "data/raw/{sample}_R1.fastq",
        r2 = "data/raw/{sample}_R2.fastq"
    output:
        r1 = "data/trimmed/{sample}_R1.fastq",
        r2 = "data/trimmed/{sample}_R2.fastq"
    shell:
        """
        trimmomatic PE {input.r1} {input.r2} {output.r1} {output.r2}
        """
```

`{sample}` - это wildcard! Snakemake сам понимает, какие значения туда подставить.

## Как это работает?

1. Snakemake смотрит на `rule all` - какие файлы нужны на выходе
2. Смотрит на существующие файлы в директориях
3. Сопоставляет паттерны и понимает значения wildcards
4. Запускает правила для каждого уникального значения

## Пример с несколькими wildcards

```python
rule align:
    input:
        reads = "data/trimmed/{sample}_{rep}_R{read}.fastq",
        reference = "data/reference.fasta"
    output:
        bam = "data/aligned/{sample}_{rep}.bam"
    shell:
        """
        bwa mem {input.reference} {input.reads} > {output.bam}
        """
```

Здесь три wildcards: `{sample}`, `{rep}`, `{read}`

## Микро-задачка 🎯

**Дано:** 
- Файлы: `SRR001_R1.fastq`, `SRR001_R2.fastq`, `SRR002_R1.fastq`, `SRR002_R2.fastq`
- Нужно запустить FastQC для каждого файла

**Напишите правило с wildcards, которое:**
1. Принимает любой файл вида `{sample}_R{read}.fastq`
2. Создает отчет `{sample}_R{read}_fastqc.html`

```python
rule all:
    input:
        expand("results/{sample}_R{read}_fastqc.html", 
               sample=["SRR001", "SRR002"], 
               read=[1, 2])

rule fastqc:
    input:
        # ??? ВАШ КОД
    output:
        # ??? ВАШ КОД
    shell:
        """
        fastqc {input} -o results/
        """
```

## Функция expand()

`expand()` - генерирует все комбинации wildcards

```python
# Вместо:
["data/sample1.txt", "data/sample2.txt", "data/sample3.txt"]

# Пишем:
expand("data/{sample}.txt", sample=["sample1", "sample2", "sample3"])
```

**Несколько wildcards:**
```python
expand("data/{sample}_R{read}.fastq", 
       sample=["SRR001", "SRR002"], 
       read=[1, 2])
# Получим:
# ["data/SRR001_R1.fastq", "data/SRR001_R2.fastq",
#  "data/SRR002_R1.fastq", "data/SRR002_R2.fastq"]
```

## Ограничение wildcards (constraints)

Иногда нужно ограничить, что может подставляться в wildcard:

```python
rule process:
    input:
        "data/{sample,[A-Z]+\d+}.txt"  # Только заглавные буквы + цифры
    output:
        "results/{sample}.txt"
    shell:
        "process.py {input} > {output}"
```

Регулярные выражения в wildcards помогают избежать конфликтов!