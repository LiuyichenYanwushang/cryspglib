#!/usr/bin/env python3
import csv, json, sys
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.iso_irrep_exact import load_exact_iso_irrep_sources
from scripts import generate_irrep_data as g
from scripts import check_other_wave_vector_rows as w

def rank(vectors):
    a=[[Fraction(x) for x in v] for v in vectors if v is not None and any(v)]
    p=0
    for c in range(3):
        q=next((i for i in range(p,len(a)) if a[i][c]),None)
        if q is None: continue
        a[p],a[q]=a[q],a[p]
        x=a[p][c]; a[p]=[v/x for v in a[p]]
        for i in range(len(a)):
            if i!=p:
                x=a[i][c]; a[i]=[v-x*b for v,b in zip(a[i],a[p])]
        p+=1
    return p

def main():
    base=Path(sys.argv[1] if len(sys.argv)>1 else 'target/irrep')
    db=load_exact_iso_irrep_sources()
    summary={}
    for name, records in [('pir',db.pir_records),('cir',db.cir_records)]:
        counts=Counter(); witnesses=defaultdict(list)
        with open(str(base)+'-'+name+'-families.tsv','w') as f:
            writer=csv.writer(f,delimiter='\t'); writer.writerow(['irnumber','sg','label','dimension','rank','constant','parameters','stored_arms','generic_kcount','irtype','operations','irtranslations'])
            for r in records:
                k=r.k_arms[0]; d=rank(k.parameters)
                assert all(rank(a.parameters)==d for a in r.k_arms)
                key=str(d); counts[key]+=1
                if d==0: counts['gamma' if not any(k.constant) else 'fixed_non_gamma']+=1
                if d and any(k.constant): counts['rank_'+key+'_nonzero_offset']+=1
                row=[r.irnumber,r.spacegroup,r.irrep_label,r.dimension,d,str(k.constant),str(k.parameters),len(r.k_arms),r.kcount,r.irtype,len(r.operations),sum(t is not None for t in r.irtranslations)]
                writer.writerow(row)
                if len(witnesses[key])<8: witnesses[key].append(row)
        summary[name]={'total':len(records),'counts':dict(counts),'witnesses':dict(witnesses)}

    # Cross-check runtime ordinary row source ids against authoritative CIR affine rank.
    rows=list(csv.DictReader(open(str(base)+'-inventory.tsv'),delimiter='\t'))
    rr=Counter(); examples=defaultdict(list)
    for row in rows:
        key=row['family']
        rr[key]+=1;rr[key+'/point_'+row['is_point']]+=1
        rr[key+'/matrix_'+str(int(row['matrix_count'])>0)]+=1
        if row['family']!='spinor':
            ranks=[]
            for irn in map(int,row['cir_sources'].split(',')):
                source=db.cir_records[irn-1]; assert source.irnumber==irn
                ranks.append(rank(source.k_arms[0].parameters))
            assert len(set(ranks))==1
            d=ranks[0]; rr[key+'/affine_rank_'+str(d)]+=1
            if d and len(examples['runtime_parameterized'])<8: examples['runtime_parameterized'].append(row)
    summary['runtime']=dict(rr)
    summary['runtime_examples']=dict(examples)

    g._verify_pinned_archives()
    data=w.parse_pinned(); resolved,failures=w.little_wave_vectors(data,list(zip(data['labels'],data['source_sg'])))
    assert not failures, failures
    summary['frozen_w']={'sources':len(resolved),'sgs':sorted(set(x[0] for x in resolved)),'parameter_counts':dict(Counter(len(x[-1]) for x in resolved)),'basepoints':dict(Counter(x[-2] for x in resolved)),'records':sum(c>0 for c in data['count']),'rows':sum(data['count'])}
    # Affine domains in official little table. Count actual occupied slots, not buffer capacity.
    lines=g.read_file('data_little.txt'); sec=g.get_sections(lines)
    kdims=g.parse_ints(lines,sec,'little_k_dim')
    assert len(kdims)==14*27
    little=Counter(); little_witnesses=defaultdict(list)
    for lattice,n in enumerate(data['little_k_count']):
        for slot in range(n):
            groups=w.wave_vector(data['little_k'],lattice,slot)
            vectors=[tuple(Fraction(x,vec[3]) for x in vec[:3]) if vec[3] else None for vec in groups[1:]]
            d=rank(vectors); assert d==kdims[lattice*27+slot],(lattice,slot,d,kdims[lattice*27+slot])
            little[d]+=1
            if len(little_witnesses[d])<4: little_witnesses[d].append({'lattice':lattice+1,'slot':slot+1,'label':data['little_k_label'][lattice*27+slot],'groups':groups})
    summary['little_domains']={'total':sum(little.values()),'ranks':dict(little),'witnesses':dict(little_witnesses)}
    # Official full little-irrep rows, joined by (SG, lattice, k-slot).
    lir=Counter(); lw=defaultdict(list)
    for sg,label,slot,dim in zip(data['little_sg'],data['little_label'],data['little_k_index'],data['little_dim']):
        if sg==0:
            assert label=='' and slot==0 and dim==0
            continue
        lattice=data['sg_lattice'][sg-1]-1; d=kdims[lattice*27+slot-1]; lir[d]+=1
        if len(lw[d])<4: lw[d].append([sg,label,dim,lattice+1,slot])
    summary['little_full_irreps']={'total':sum(lir.values()),'unused_buffer_rows':sum(sg==0 for sg in data['little_sg']),'ranks':dict(lir),'witnesses':dict(lw)}
    assert len(rows)==8388 and len({r['id'] for r in rows})==8388
    assert all(r['view']=='ok' for r in rows)
    assert all(r['dim']==r['star_dim'] for r in rows if r['spinor']=='false')
    assert {k:summary['pir']['counts'][k] for k in ['0','1','2','3']}=={'0':4777,'1':3980,'2':1307,'3':230}
    assert summary['little_full_irreps']['ranks']=={0:4777,1:3980,2:1307,3:230}
    assert all(int(r['matrix_count'])==0 for r in rows if r['spinor']=='true')
    assert len(resolved)==73 and sum(data['count'])==5756
    Path(str(base)+'-archives-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__ == "__main__":
    main()
