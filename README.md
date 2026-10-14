# bioinf-utils

Набор утилит для работы с биоинформатическими форматами: FASTQ, FASTA, BLAST.

# Установка

git clone https://github.com/DeniKrav/bioinf-utils.git
cd bioinf-utils

Требуется Python 3.8+. Сторонние библиотеки не нужны.

# Использование

# Фильтрация FASTQ

from bioinf_utils import filter_fastq

filter_fastq(
    "example_data/reads.fastq",
    "filtered.fastq",
    gc_bounds=(40, 60),
    length_bounds=(50, 150),
    quality_threshold=20,
)

Результат сохраняется в папку filtered/ под именем filtered.fastq.
gc_bounds и length_bounds можно передать одним числом — тогда это верхняя граница.
Границы включаются. quality_threshold — Phred33.

# FASTA: многострочный → односточный

from bioinf_utils import convert_multiline_fasta_to_oneline

convert_multiline_fasta_to_oneline(
    "example_data/multiline.fasta",
    "example_data/oneline.fasta",
)

Без второго аргумента имя формируется автоматически: <input>_oneline.fasta.

# BLAST: извлечь лучшие совпадения

from bioinf_utils import parse_blast_output

parse_blast_output(
    "example_data/example_blast_results.txt",
    "example_data/proteins.txt",
)

Для каждого Query= берётся первая строка из блока Sequences producing significant alignments:.
Результат — отсортированный по алфавиту список белков, по одному в строке.

# Структура

bioinf-utils/
├── bioinf_utils.py       # публичные функции
├── helpers/
│   ├── __init__.py
│   ├── fastq.py          # вспомогательные функции для FASTQ
│   └── bio_files.py      # вспомогательные функции для FASTA и BLAST
└── example_data/         # примеры входных файлов

# Функции

filter_fastq                        — фильтрация FASTQ по GC, длине, качеству
convert_multiline_fasta_to_oneline  — многострочный FASTA → односточный
parse_blast_output                  — извлечение лучших совпадений из BLAST txt

# Контакты

Денис Кравецкий
GitHub: https://github.com/DeniKrav