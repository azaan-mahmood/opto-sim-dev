"""Generate LaTeX table .tex files from validation _table.csv files.

Usage:
    python generate_latex_tables.py [--seed 42]

Reads *_table.csv from each val_*/ output directory and writes
corresponding .tex files to paperwork/tables/.
"""
import csv
import os
import argparse

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(SCRIPT_DIR, '..', '..')
TABLES_OUT = os.path.join(BASE, 'paperwork', 'tables')
os.makedirs(TABLES_OUT, exist_ok=True)

# Map: (csv_subdir, csv_prefix, tex_label, caption)
TABLES = [
    ('val_cwlaser',       'val_cwlaser',       'tab:cwlaser',       'CW laser validation summary.'),
    ('val_mzm',           'val_mzm',            'tab:mzm',           'MZM validation summary.'),
    ('val_cd',            'val_cd',             'tab:cd',            'Chromatic dispersion validation summary.'),
    ('val_pmd',           'val_pmd',            'tab:pmd',           'PMD validation summary.'),
    ('val_attenuation',   'val_attenuation',    'tab:attenuation',   'Attenuation validation summary.'),
    ('val_birefringence', 'val_birefringence',  'tab:birefringence', 'Birefringence validation summary.'),
    ('val_apd',           'val_apd',            'tab:apd',           'APD validation summary.'),
    ('val_system',        'val_system',         'tab:system',        'System-level combined impairment scenarios at 100~km.'),
]


def escape_latex(s):
    """Escape special LaTeX characters in a string."""
    s = s.replace('\\', '\\textbackslash{}')
    for ch in '&%$#_{}':
        s = s.replace(ch, '\\' + ch)
    s = s.replace('~', '\\textasciitilde{}')
    s = s.replace('^', '\\textasciicircum{}')
    return s


def read_csv(csv_path):
    """Read a CSV file and return (headers, rows)."""
    with open(csv_path, 'r') as f:
        reader = csv.reader(f)
        headers = next(reader)
        rows = list(reader)
    return headers, rows


def guess_col_spec(n_cols):
    """Return a tabular column spec for n_cols columns."""
    if n_cols == 2:
        return 'l l'
    elif n_cols == 3:
        return 'l c c'
    elif n_cols == 4:
        return 'c c c c'
    return 'l ' + 'c ' * (n_cols - 1)


def generate_tex(headers, rows, label, caption, n_cols):
    """Generate a LaTeX table environment string."""
    col_spec = guess_col_spec(n_cols)
    lines = [
        '\\begin{table}[!ht]',
        '  \\caption{' + caption + '}',
        '  \\label{' + label + '}',
        '  \\centering',
        '  \\begin{tabular}{' + col_spec + '}',
        '    \\toprule',
        '    ' + ' & '.join(escape_latex(h) for h in headers) + ' \\\\',
        '    \\midrule',
    ]
    for row in rows:
        lines.append('    ' + ' & '.join(escape_latex(cell) for cell in row) + ' \\\\')
    lines.append('    \\bottomrule')
    lines.append('  \\end{tabular}')
    lines.append('\\end{table}')
    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()
    seed = args.seed

    for subdir, csv_prefix, tex_label, caption in TABLES:
        csv_path = os.path.join(BASE, 'analysis', subdir if subdir != 'val_system' else '',
                                f'{csv_prefix}--seed{seed}_table.csv')
        if not os.path.exists(csv_path):
            # Try alternate path for val_system
            csv_path = os.path.join(BASE, subdir, f'{csv_prefix}--seed{seed}_table.csv')
        if not os.path.exists(csv_path):
            print(f"  SKIP: {csv_path} not found")
            continue

        headers, rows = read_csv(csv_path)
        n_cols = len(headers)
        tex_content = generate_tex(headers, rows, tex_label, caption, n_cols)

        tex_path = os.path.join(TABLES_OUT, f'{csv_prefix}_table.tex')
        with open(tex_path, 'w') as f:
            f.write(tex_content)
        print(f"  Wrote: {tex_path}")

    print("Done.")


if __name__ == '__main__':
    main()
