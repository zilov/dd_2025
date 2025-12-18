# Python argparse

## Зачем нужен argparse?

**Проблема:** Хардкод параметров в скрипте - плохая практика
```python
input_file = "data.txt"  # А если нужен другой файл?
threads = 8              # А если на сервере другое количество?
```

**Решение:** argparse - делает ваш скрипт гибким и профессиональным

## Базовый пример

```python
import argparse

parser = argparse.ArgumentParser(description='Process sequencing data')

parser.add_argument('--input', '-i', 
                    required=True,
                    help='Input FASTQ file')

parser.add_argument('--output', '-o',
                    default='output.txt',
                    help='Output file')

parser.add_argument('--threads', '-t',
                    type=int,
                    default=4,
                    help='Number of threads')

args = parser.parse_args()

# Теперь используем:
print(f"Processing {args.input}")
print(f"Using {args.threads} threads")
```

**Запуск:**
```bash
python script.py --input data.fastq --threads 8
python script.py -i data.fastq -t 8  # короткая форма
```

## Типы аргументов

### Позиционные (обязательные)
```python
parser.add_argument('input_file')  # Без -- это позиционный аргумент
```
Запуск: `python script.py myfile.txt`

### Опциональные
```python
parser.add_argument('--output', default='out.txt')
```

### Флаги (boolean)
```python
parser.add_argument('--verbose', '-v', action='store_true')
# Если указан --verbose, то args.verbose = True
```

### С выбором из вариантов
```python
parser.add_argument('--format', choices=['fastq', 'fasta', 'bam'])
```

## Типы данных

```python
parser.add_argument('--threads', type=int)
parser.add_argument('--quality', type=float)
parser.add_argument('--samples', nargs='+')  # Список значений
```

**Запуск:**
```bash
python script.py --threads 8 --quality 0.95 --samples S1 S2 S3
```

## Автоматическая справка

argparse автоматически создает `-h` и `--help`:

```bash
python script.py --help
```

Выведет:
```
usage: script.py [-h] --input INPUT [--output OUTPUT] [--threads THREADS]

Process sequencing data

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        Input FASTQ file
  --output OUTPUT, -o OUTPUT
                        Output file
  --threads THREADS, -t THREADS
                        Number of threads
```

## Практический пример для биоинформатики

```python
import argparse

def main():
    parser = argparse.ArgumentParser(
        description='Run FastQC quality control'
    )
    
    parser.add_argument('--reads', '-r',
                        nargs='+',
                        required=True,
                        help='Input FASTQ files')
    
    parser.add_argument('--outdir', '-o',
                        default='fastqc_results',
                        help='Output directory')
    
    parser.add_argument('--threads', '-t',
                        type=int,
                        default=1,
                        help='Number of threads')
    
    parser.add_argument('--quiet', '-q',
                        action='store_true',
                        help='Suppress output')
    
    args = parser.parse_args()
    
    # Ваш код здесь
    for read_file in args.reads:
        print(f"Processing {read_file}")
        # run_fastqc(read_file, args.outdir, args.threads)

if __name__ == '__main__':
    main()
```

## Интеграция со Snakemake

```python
# В Snakefile можно вызвать ваш скрипт:
rule run_analysis:
    input:
        reads="data/{sample}.fastq"
    output:
        result="results/{sample}.txt"
    threads: 8
    shell:
        """
        python scripts/analyze.py \
            --input {input.reads} \
            --output {output.result} \
            --threads {threads}
        """
```

## Best Practices

1. **Всегда добавляйте help** - коллеги скажут спасибо
2. **Используйте значения по умолчанию** - для частых случаев
3. **Валидируйте входные данные** - проверяйте существование файлов
4. **Группируйте аргументы** - для сложных программ используйте `add_argument_group()`

```python
parser = argparse.ArgumentParser()
input_group = parser.add_argument_group('Input options')
input_group.add_argument('--reads', required=True)
input_group.add_argument('--reference')

output_group = parser.add_argument_group('Output options')
output_group.add_argument('--outdir', default='results')
```