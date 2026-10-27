#!/usr/bin/env python3
"""Audit the frozen all-UNI census; no production code is imported or changed."""
import csv, json, collections, pathlib, hashlib, re, math, sys
P=pathlib.Path(sys.argv[1] if len(sys.argv)>1 else 'target/magnetic-target-')
def path(s): return pathlib.Path(str(P)+s)
def load(s):
    with path(s+'.tsv').open() as f: return list(csv.DictReader(f,delimiter='\t'))
u,g,s,x=(load(n) for n in ('uni','gamma','sources','unresolved'))
assert len(u)==1651 and {int(r['uni']) for r in u}==set(range(1,1652))
assert len(g)==sum(int(r['gamma_targets']) for r in u)
assert len(x)==sum(int(r['gamma_unresolved'])+int(r['other_unresolved']) for r in u)
assert len({(r['uni'],r['index']) for r in g})==len(g)
assert all(r['status']=='OK' and r['gamma_points']=='1' for r in u)
sm={r['key']:r for r in s}; um={r['uni']:r for r in u}
assert len(sm)==len(s)
assert all(r['named_owners']=='1' and r['authoritative_rows']=='1' for r in s)
assert all('ERROR=' not in r['detail'] and 'NO_EXPOSED_MATCH' not in r['detail'] for r in s)
owner_ids={re.search(r'owner=.*?#(\d+)\[',r['detail']).group(1) for r in s}
ordinary_owners={re.search(r'owner=.*?#(\d+)\[',r['detail']).group(1) for r in s if r['spinor']=='false'}
assert all(('raw_elements=0,' in r['detail']) if r['spinor']=='true' else ('raw_shape=true' in r['detail'] and 'finite=true' in r['detail'] and 'typed_trace=Some(true)' in r['detail']) for r in s)
for r in g:
    unitary=int(r['anti_order'])==0
    assert unitary==r['family'].startswith('unitary-')
    assert r['raw_type'] in ('A','B','C')
    if unitary:
        assert r['raw_type']=='A'
        assert r['family']==('unitary-packet' if int(r['h_source_count'])>1 else 'unitary-single')
    else:
        assert r['family']=='magnetic-'+r['raw_type']
        assert int(r['h_source_count'])==(2 if r['raw_type']=='C' else 1)
    expected_norm=4 if r['raw_type']=='B' else 2 if r['raw_type']=='C' or r['family']=='unitary-packet' else 1
    assert math.isclose(float(r['unitary_norm']),expected_norm,abs_tol=1e-8)
    hs=[sm[k] for k in r['h_keys'].split(';')]
    assert len(hs)==int(r['h_source_count'])
    assert all(h['spinor']==r['spinor'] and h['h_sg']==r['h_sg'] for h in hs)
    r['h_storage']='+'.join(h['status'] for h in hs)
    r['h_complex']='+'.join('complex' if 'complex=true' in h['detail'] else 'real' for h in hs)
    r['h_available']=str(all(h['authoritative_rows']=='1' for h in hs))
    r['isotropy_record_count']=um[r['uni']]['isotropy_records']
    r['record_witness']=um[r['uni']]['isotropy_witness']
    r['scope_family']=('spinor' if r['spinor']=='true' else 'nonspinor-B' if r['family']=='magnetic-B' else 'nonspinor-3D-A' if r['dimension']=='3' else 'existing-target-family-shape')
def counter(rows,keys):return dict(sorted(collections.Counter('|'.join(r[k] for k in keys) for r in rows).items()))
def firsts(rows,keys):
    out={}
    for r in rows:out.setdefault('|'.join(r[k] for k in keys),r)
    return out
counts={
 'denominators':{'uni':len(u),'gamma_summaries':len(g),'other_k_summaries':sum(int(r['other_targets']) for r in u),'unique_H_sources':len(s),'H_source_occurrences':sum(int(r['h_source_count']) for r in g),'unresolved_Gamma':sum(int(r['gamma_unresolved']) for r in u),'unresolved_other_k':sum(int(r['other_unresolved']) for r in u),'fatal':sum(r['status']!='OK' for r in u),'isotropy_records':sum(int(r['isotropy_records']) for r in u),'isotropy_UNIs':sum(int(r['isotropy_records'])>0 for r in u)},
 'distinct_physical_H_owners':{'all':len(owner_ids),'nonspinor':len(ordinary_owners),'spinor':len(owner_ids-ordinary_owners)},
 'uni_type':counter(u,['magnetic_type']),
 'uni_missing_isotropy':counter([r for r in u if r['isotropy_records']=='0'],['magnetic_type']),
 'Gamma_raw_type':counter(g,['raw_type','spinor']),
 'Gamma_disjoint_family':counter(g,['family','spinor','dimension','h_source_count','h_dimensions']),
 'Gamma_storage':counter(g,['family','spinor','dimension','h_storage','h_complex']),
 'Gamma_scope_shape':counter(g,['scope_family']),
 'Gamma_completeness':counter([{**r,'completion':r['completeness'].split(' ')[0]} for r in g],['completion','family','spinor']),
 'Gamma_empty_candidate_lists':counter([r for r in g if r['ordinary_candidates']=='0' and r['magnetic_candidates']=='0'],['family','spinor','h_storage']),
 'Gamma_shape_witnesses':firsts(g,['scope_family']),
 'unique_H_storage':counter(s,['status','dimension']),
 'unique_H_owner_cardinality':counter(s,['named_owners','authoritative_rows','standalone_rows']),
 'H_source_occurrence_storage':counter([sm[k] for r in g for k in r['h_keys'].split(';')],['status','dimension']),
 'B_nonspinor_UNIs':sorted({int(r['uni']) for r in g if r['family']=='magnetic-B' and r['spinor']=='false'}),
 'unresolved':x,
 'sha256':{n:hashlib.sha256(path(n+'.tsv').read_bytes()).hexdigest() for n in ('uni','gamma','sources','unresolved')},
}
assert len(g)==12627 and len(s)==2128 and len(x)==0
assert counts['Gamma_scope_shape']=={'existing-target-family-shape':8037,'nonspinor-3D-A':320,'nonspinor-B':20,'spinor':4250}
assert counts['distinct_physical_H_owners']=={'all':2059,'nonspinor':1313,'spinor':746}
assert counts['denominators']['other_k_summaries']==40166
assert counts['denominators']['isotropy_records']==16721 and counts['denominators']['isotropy_UNIs']==1421
assert counts['uni_missing_isotropy']=={'Grey':230}
path('counts.json').write_text(json.dumps(counts,indent=2)+'\n')
with path('classified.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(g[0]),delimiter='\t');w.writeheader();w.writerows(g)
print(json.dumps(counts,indent=2))
