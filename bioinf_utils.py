"""Биоинформатические утилиты: FASTQ-фильтрация, работа с GBK и BLAST."""

import os

from helpers.fastq import (
    read_fastq,
    write_fastq_record,
    passes_gc,
    passes_length,
    passes_quality,
)
from helpers.bio_files import (
    read_fasta,
    write_fasta_record,
    read_blast_txt,
    extract_first_description,
    write_lines,
)


def filter_fastq(input_fastq, output_fastq,
                 gc_bounds=(0, 100),
                 length_bounds=(0, 2**32),
                 quality_threshold=0):
    """Фильтрует FASTQ-файл по GC-составу, длине и среднему качеству.

    Отфильтрованные риды сохраняются в папку filtered/ под именем output_fastq.
    Работает «на лету»: в памяти — не более одной записи.

    Args:
        input_fastq: путь до входного fastq-файла.
        output_fastq: имя выходного файла (сохраняется в filtered/).
        gc_bounds: интервал GC-состава в процентах (0..100).
            Одно число — верхняя граница.
        length_bounds: интервал длины рида. Одно число — верхняя граница.
        quality_threshold: минимальное среднее качество (Phred33).
    """
    os.makedirs("filtered", exist_ok=True)
    output_path = os.path.join("filtered", output_fastq)

    with open(input_fastq) as fin, open(output_path, "w") as fout:
        for header, seq, qual in read_fastq(fin):
            if not passes_length(seq, length_bounds):
                continue
            if not passes_gc(seq, gc_bounds):
                continue
            if not passes_quality(qual, quality_threshold):
                continue
            write_fastq_record(fout, header, seq, qual)


def convert_multiline_fasta_to_oneline(input_fasta, output_fasta=None):
    """Преобразует многострочный FASTA в односточный.

    Читает FASTA-файл, где последовательность может быть разбита
    на несколько строк, и сохраняет новый FASTA, где каждая
    последовательность записана одной строкой.

    Args:
        input_fasta: путь до входного FASTA-файла.
        output_fasta: имя выходного файла. Если не задано — формируется
            автоматически как <input>_oneline.fasta.
    """
    if output_fasta is None:
        base, _ = os.path.splitext(input_fasta)
        output_fasta = base + "_oneline.fasta"

    with open(input_fasta) as fin, open(output_fasta, "w") as fout:
        for header, sequence in read_fasta(fin):
            write_fasta_record(fout, header, sequence)


def parse_blast_output(input_file, output_file):
    """Извлекает первые Description для каждого QUERY из BLAST-файла.

    Читает txt-файл с результатами BLAST, для каждого запроса берёт
    первую строку из блока 'Sequences producing significant alignments:'
    и сохраняет список полученных белков в output_file, отсортированный
    по алфавиту, по одному в строке.

    Args:
        input_file: путь до входного txt-файла с результатами BLAST.
        output_file: путь до выходного файла (по одному белку в строке).
    """
    blocks = read_blast_txt(input_file)

    proteins = []
    for block in blocks:
        desc = extract_first_description(block)
        if desc is not None:
            proteins.append(desc)

    proteins.sort()
    write_lines(output_file, proteins)