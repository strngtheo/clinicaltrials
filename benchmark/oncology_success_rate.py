import csv
import ast
from collections import defaultdict
from ccs_utils import cancer_filter_icd10code


def parse_icd_codes(text):
    """Parse the ICD codes column into a list of codes without dots."""
    try:
        outer = ast.literal_eval(text)
    except Exception:
        return []
    codes = []
    for item in outer:
        try:
            inner = ast.literal_eval(item) if isinstance(item, str) else item
            for code in inner:
                codes.append(code.replace('.', ''))
        except Exception:
            continue
    return codes


def parse_diseases(text):
    try:
        return [d.lower() for d in ast.literal_eval(text)]
    except Exception:
        return []


def is_oncology(disease_text, icd_text):
    diseases = parse_diseases(disease_text)
    if any(keyword in ' '.join(diseases) for keyword in ['cancer', 'neoplasm', 'tumor']):
        return True
    for code in parse_icd_codes(icd_text):
        try:
            if cancer_filter_icd10code(code):
                return True
        except KeyError:
            continue
    return False


def load_start_years(path):
    """Return a mapping from nctid to start year."""
    mapping = {}
    with open(path) as f:
        for line in f:
            nctid, start_date, _ = line.strip().split('\t')
            year = 0
            if start_date:
                try:
                    year = int(start_date.split()[-1])
                except ValueError:
                    year = 0
            mapping[nctid] = year
    return mapping


def oncology_success_rates(raw_csv='data/raw_data.csv', date_file='data/nctid_date.txt'):
    years = load_start_years(date_file)
    stats = defaultdict(lambda: [0, 0])
    with open(raw_csv) as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            nctid = row[0]
            label = int(row[3])
            if nctid not in years or years[nctid] == 0:
                continue
            if is_oncology(row[5], row[6]):
                y = years[nctid]
                stats[y][0] += label
                stats[y][1] += 1
    rates = {y: stats[y][0] / stats[y][1] * 100 for y in stats if stats[y][1] > 0}
    return rates


def main():
    rates = oncology_success_rates()
    for year in sorted(rates):
        print(f"{year}\t{rates[year]:.2f}")


if __name__ == "__main__":
    main()
