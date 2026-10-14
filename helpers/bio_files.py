"""Вспомогательные функции для работы с биоинформатическими файлами."""


def read_fasta(path):
    """Генератор: читает FASTA по одной записи.

    Последовательность может быть разбита на несколько строк —
    они склеиваются в одну.

    Yields:
        (header, sequence) — заголовок (со знаком '>') и склеенная последовательность.
    """
    header = None
    parts = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line
                parts = []
            else:
                parts.append(line.strip())
    if header is not None:
        yield header, "".join(parts)


def write_fasta_record(out, header, sequence):
    """Записывает одну FASTA-запись в открытый файл.

    Args:
        out: открытый на запись файловый объект.
        header: строка заголовка (со знаком '>').
        sequence: строка последовательности (без переводов строк).
    """
    out.write(f"{header}\n{sequence}\n")


def read_blast_txt(path):
    """Читает BLAST txt, возвращает список блоков — по одному на QUERY.

    Каждый блок — список строк между 'Query=' и следующим 'Query='
    (или концом файла).
    """
    blocks = []
    current = None
    with open(path) as f:
        for line in f:
            if line.startswith("Query="):
                if current is not None:
                    blocks.append(current)
                current = []
            elif current is not None:
                current.append(line)
    if current is not None:
        blocks.append(current)
    return blocks


def extract_first_description(block):
    """Из блока BLAST вытаскивает Description первой строки совпадений.

    Возвращает строку (название белка) или None, если блок пустой.
    """
    in_alignments = False
    header_skipped = False

    for line in block:
        if line.startswith("Sequences producing significant alignments:"):
            in_alignments = True
            continue
        if not in_alignments:
            continue

        stripped = line.strip()
        if not stripped:
            continue

        if not header_skipped:
            header_skipped = True
            continue

        return extract_description(stripped)

    return None


def extract_description(line):
    """Из строки BLAST-хита извлекает Description (первый столбец).

    Столбцы разделены двумя и более пробелами. Description может
    содержать одинарные пробелы, поэтому split("  ") безопасен.
    """
    parts = line.split("  ")
    parts = [p.strip() for p in parts if p.strip()]
    return parts[0] if parts else None


def write_lines(path, lines):
    """Записывает список строк в файл — по одной в строке."""
    with open(path, "w") as f:
        for line in lines:
            f.write(line + "\n")