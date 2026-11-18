#!/usr/bin/env python3
"""Generate reviewable metadata from chromosome-arm .seq filenames.

Default: the leading alphabetical prefix immediately preceding the first digit
is the population (CK1 -> CK; EA49 -> EA; ZI90 -> ZI). Unknown schemes require
an explicit mapping, never a guessed population.
"""
import argparse
import csv
import re
from pathlib import Path

ARMS = r'2L|2R|3L|3R|X'
FILE_RE = re.compile(rf'(?P<sample>.+)_Chr(?P<arm>{ARMS})\.seq$', re.I)


def read_mapping(path, columns):
    if not path: return {}
    with open(path, newline='') as handle:
        r = csv.DictReader(handle)
        if not r.fieldnames or not set(columns) <= set(r.fieldnames):
            raise ValueError(f'{path} requires columns {",".join(columns)}')
        out = {}
        for row in r:
            key = row[columns[0]].strip()
            val = row[columns[1]].strip()
            if not key or not val: raise ValueError(f'Empty mapping: {row}')
            if key in out and out[key] != val: raise ValueError(f'Conflicting mapping for {key}')
            out[key] = val
        return out


def infer_prefix(sample, expression):
    match = re.match(expression, sample)
    if not match:
        return None
    return match.group('population') if 'population' in match.groupdict() else match.group(1)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--seq-dir', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--sample-map', type=Path, help='CSV sample,population; overrides inferred labels')
    p.add_argument('--prefix-map', type=Path, help='CSV prefix,population; overrides inferred prefix mapping')
    p.add_argument('--prefix-regex', default=r'^([A-Za-z]+)(?=\d)',
                   help='Regex to infer population from sample name; first group or named group population')
    p.add_argument('--case', choices=['upper','lower','preserve'], default='upper',
                   help='Case normalization of inferred population labels, not explicit mappings')
    p.add_argument('--unknown', choices=['error','skip'], default='error',
                   help='How to handle ambiguous sample naming without mappings')
    args = p.parse_args(argv)
    by_sample = {}
    for f in args.seq_dir.iterdir():
        match = FILE_RE.fullmatch(f.name)
        if f.is_file() and match:
            by_sample.setdefault(match.group('sample'), set()).add(match.group('arm').upper())
    if not by_sample: p.error('No matching *_Chr2L.seq / *_Chr3R.seq files found')
    smap = read_mapping(args.sample_map, ('sample','population'))
    pmap = read_mapping(args.prefix_map, ('prefix','population'))
    pmap = {k.upper():v for k,v in pmap.items()}
    missing = []
    rows=[]
    for sample in sorted(by_sample):
        prefix = infer_prefix(sample,args.prefix_regex)
        if sample in smap:
            pop = smap[sample]
        elif prefix and prefix.upper() in pmap:
            pop = pmap[prefix.upper()]
        elif prefix:
            pop = {'upper':prefix.upper(),'lower':prefix.lower(),'preserve':prefix}[args.case]
        else:
            missing.append(sample)
            continue
        rows.append((sample,pop))
    if missing and args.unknown=='error':
        p.error(f'Could not infer population for {len(missing)} samples: {missing[:10]}. '
                'Supply --sample-map / --prefix-map or a suitable --prefix-regex. No output written.')
    if not rows: p.error('No samples assigned to populations')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w',newline='') as handle:
        writer=csv.writer(handle)
        writer.writerow(['sample','population'])
        writer.writerows(rows)
    print(f'Wrote {len(rows)} samples ({len(set(pop for _,pop in rows))} populations) to {args.output}')
    if missing:print(f'Skipped {len(missing)} samples without assignable population: {missing[:10]}')
    print('Population sizes:')
    from collections import Counter
    for pop,n in sorted(Counter(pop for _,pop in rows).items()):print(f'  {pop}: {n}')
    print('Review assignments before interpreting results. Arms per sample:')
    for sample,pop in rows:print(f'{sample}\t{pop}\t{",".join(sorted(by_sample[sample]))}')


if __name__ == '__main__':main()
