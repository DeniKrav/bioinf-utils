"""Вспомогательные функции для работы с FASTQ."""


def read_fastq(path):
    """Генератор: читает FASTQ по одной записи.

    Yields:
        (header, seq, qual) — три строки без переводов строк.
    """
    with open(path) as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()          # строка "+"
            qual = f.readline().strip()
            yield header, seq, qual


def write_fastq_record(out, header, seq, qual):
    """Записывает одну FASTQ-запись в открытый файл out."""
    out.write(f"{header}\n{seq}\n+\n{qual}\n")


def gc_content(seq):
    """GC-состав в процентах (0..100)."""
    if not seq:
        return 0.0
    gc = seq.upper().count("G") + seq.upper().count("C")
    return gc / len(seq) * 100


def mean_quality(qual):
    """Среднее качество рида по шкале Phred33."""
    if not qual:
        return 0.0
    return sum(ord(ch) - 33 for ch in qual) / len(qual)


def normalize_bounds(bounds, default_max):
    """Приводит границы к кортежу (low, high).

    Если передано одно число X — считаем его ВЕРХНЕЙ границей: (0, X).
    Если кортеж — возвращаем как есть.
    """
    if isinstance(bounds, (int, float)):
        return (0, bounds)
    return bounds


def passes_gc(seq, gc_bounds):
    """True, если GC-состав рида попадает в интервал (включительно)."""
    low, high = normalize_bounds(gc_bounds, 100)
    gc = gc_content(seq)
    return low <= gc <= high


def passes_length(seq, length_bounds):
    """True, если длина рида попадает в интервал (включительно)."""
    low, high = normalize_bounds(length_bounds, 2**32)
    return low <= len(seq) <= high


def passes_quality(qual, quality_threshold):
    """True, если среднее качество >= порога."""
    return mean_quality(qual) >= quality_threshold