#!/usr/bin/env python3
"""Convert FlyBase GO GAF 2.x to GO term GMT and optional exploratory ion-focused GMT.

No third-party Python dependencies. Gene IDs remain FBgn stable identifiers.
Ontology propagation follows is_a and (optionally) part_of only, not regulates.
"""
import argparse
import csv
import gzip
import json
import re
from collections import defaultdict, Counter
from functools import lru_cache
from pathlib import Path

ROOTS = {'GO:0008150', 'GO:0003674', 'GO:0005575'}
# Narrow term-label patterns; they define candidates, NOT a manually curated pathway database.
ION_PATTERNS = {
    'Manganese_related': r'\bmanganese\b|\bmanganous\b',
    'Iron_related': r'\biron\b|\bferric\b|\bferrous\b',
    'Zinc_related': r'\bzinc\b',
    'Calcium_related': r'\bcalcium\b',
    'Magnesium_related': r'\bmagnesium\b',
    'Sodium_related': r'\bsodium\b',
    'Potassium_related': r'\bpotassium\b',
    'Phosphate_related': r'\bphosphate\b|\bphosphorus\b',
    'Sulfur_related': r'\bsulfur\b|\bsulphur\b|\bsulfate\b|\bsulphate\b',
    'Strontium_related': r'\bstrontium\b',
    'Metal_ion_homeostasis': r'\bmetal ion homeostasis\b|\btransition metal ion homeostasis\b',
    'Oxidative_stress': r'\boxidative stress\b|\bresponse to reactive oxygen species\b',
}


def open_text(path):
    return gzip.open(path, 'rt', encoding='utf-8') if str(path).endswith('.gz') else open(path, encoding='utf-8')


def load_obo(path, include_part_of=True):
    names, parents, altids, obsolete = {}, defaultdict(set), {}, set()
    curr = None
    def commit(t):
        if not t or 'id' not in t: return
        ident = t['id']
        if t.get('obsolete'): obsolete.add(ident)
        names[ident] = t.get('name', ident)
        parents[ident].update(t.get('parents', set()))
        for alt in t.get('alt', []): altids[alt] = ident
    with open_text(path) as inp:
        for raw in inp:
            line = raw.strip()
            if line == '[Term]':
                commit(curr); curr = {'parents': set(), 'alt': []}
            elif line.startswith('['):
                commit(curr); curr = None
            elif curr is not None:
                if line.startswith('id: '): curr['id'] = line[4:]
                elif line.startswith('name: '): curr['name'] = line[6:]
                elif line.startswith('alt_id: '): curr['alt'].append(line[8:])
                elif line.startswith('is_obsolete: true'): curr['obsolete'] = True
                elif line.startswith('is_a: '): curr['parents'].add(line[6:].split()[0])
                elif include_part_of and line.startswith('relationship: part_of '):
                    curr['parents'].add(line.split()[2])
        commit(curr)
    return names, parents, altids, obsolete


def read_gaf(path, evidence, exclude_iea, allow_db):
    direct = defaultdict(set)
    stats = Counter()
    with open_text(path) as inp:
        for line in inp:
            if not line.strip() or line.startswith('!') or line.startswith('#'): continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) < 15:
                stats['bad_columns'] += 1; continue
            db, gene, qualifier, go, ev = fields[0], fields[1], fields[3], fields[4], fields[6]
            stats['total'] += 1
            if allow_db and db not in allow_db: stats['other_db'] += 1; continue
            if 'NOT' in re.split(r'[|,]', qualifier): stats['NOT'] += 1; continue
            if exclude_iea and ev == 'IEA': stats['IEA_filtered'] += 1; continue
            if evidence and ev not in evidence: stats['evidence_filtered'] += 1; continue
            if not re.fullmatch(r'GO:\d{7}', go): stats['invalid_go'] += 1; continue
            if not gene.startswith('FBgn'): stats['non_FBgn_id'] += 1; continue
            direct[go].add(gene)
            stats['accepted'] += 1
    return direct, stats


def write_gmt(path, sets, descriptions, minimum=1):
    count=0
    with open(path, 'w', encoding='utf-8') as out:
        for name in sorted(sets):
            genes=sorted(sets[name])
            if len(genes)<minimum: continue
            desc = descriptions.get(name, 'FlyBase_GAF')
            out.write('\t'.join([name, desc.replace('\t',' ')] + genes)+'\n')
            count+=1
    return count


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--gaf', required=True, help='FlyBase gene_association.fb.gz (GAF 2.x)')
    p.add_argument('--obo', required=True, help='GO go-basic.obo (for term names and ancestor propagation)')
    p.add_argument('--out-gmt', required=True, help='Output full GO term GMT')
    p.add_argument('--ion-gmt', help='Optional exploratory ion-specific GMT based on GO term names')
    p.add_argument('--min-genes', type=int, default=5, help='Minimum distinct genes per exported set (default 5)')
    p.add_argument('--evidence', help='Comma-separated evidence codes to include; default all')
    p.add_argument('--exclude-iea', action='store_true', help='Exclude electronically inferred annotations (IEA)')
    p.add_argument('--no-propagate', action='store_true', help='Only direct annotations, do not propagate GO ancestry')
    p.add_argument('--no-part-of', action='store_true', help='Propagate is_a only (default also part_of)')
    p.add_argument('--all-db', action='store_true', help='Do not restrict GAF DB column to FB')
    p.add_argument('--report', help='Write JSON summary and ion seed term lists to this path')
    args=p.parse_args()
    if args.min_genes<1: p.error('--min-genes must be at least 1')
    names, parents, altids, obsolete=load_obo(args.obo, not args.no_part_of)
    direct, stats=read_gaf(args.gaf, set(args.evidence.split(',')) if args.evidence else None,
                           args.exclude_iea, None if args.all_db else {'FB'})

    @lru_cache(None)
    def ancestors(term):
        term=altids.get(term, term)
        if term in obsolete or term not in names: return frozenset()
        result={term}
        # clean implementation
        pending=list(parents.get(term, ()))
        while pending:
            raw=pending.pop()
            node=altids.get(raw, raw)
            if node in result or node in obsolete or node not in names: continue
            result.add(node)
            pending.extend(parents.get(node, ()))
        return frozenset(result)

    term_genes=defaultdict(set)
    skipped_unknown=0
    for term, genes in direct.items():
        term=altids.get(term,term)
        terms=ancestors(term) if not args.no_propagate else ({term} if term in names and term not in obsolete else set())
        if not terms: skipped_unknown+=len(genes); continue
        for t in terms: term_genes[t].update(genes)
    for root in ROOTS: term_genes.pop(root, None)
    descriptions={term:names[term] for term in term_genes}
    Path(args.out_gmt).parent.mkdir(parents=True,exist_ok=True)
    full_n=write_gmt(args.out_gmt, term_genes, descriptions,args.min_genes)

    report={'GAF':str(args.gaf),'OBO':str(args.obo),'stats':dict(stats),
            'unknown_obsolete_annotations':skipped_unknown,
            'all_gene_sets_exported':full_n,
            'all_terms_with_annotations':len(term_genes),
            'distinct_genes':len(set().union(*term_genes.values())) if term_genes else 0,
            'propagate':not args.no_propagate,'part_of':not args.no_part_of,
            'min_genes':args.min_genes,'ion_seed_terms':{}}
    if args.ion_gmt:
        ion_sets, ion_desc={},{}
        for label, pattern in ION_PATTERNS.items():
            rex=re.compile(pattern, re.I)
            seeds=[term for term,name in names.items() if term not in obsolete and rex.search(name)]
            genes=set().union(*(term_genes.get(t,set()) for t in seeds)) if seeds else set()
            ion_sets[label]=genes
            ion_desc[label]='GO_label_keyword_derived_REVIEW_BEFORE_USE'
            report['ion_seed_terms'][label]=[{'go_id':t,'name':names[t]} for t in sorted(seeds)]
        Path(args.ion_gmt).parent.mkdir(parents=True,exist_ok=True)
        report['ion_sets_exported']=write_gmt(args.ion_gmt,ion_sets,ion_desc,args.min_genes)
    if args.report:
        with open(args.report,'w') as out: json.dump(report,out,indent=2)
    print(f"GAF rows accepted: {stats['accepted']:,}; unique FBgn genes: {report['distinct_genes']:,}")
    print(f"GO GMT sets written: {full_n:,} -> {args.out_gmt}")
    if args.ion_gmt:print(f"Ion-focused sets written: {report['ion_sets_exported']:,} -> {args.ion_gmt}")
    print('IMPORTANT: review ion keyword-derived sets before treating them as curated biological hypotheses.')

if __name__=='__main__': main()
