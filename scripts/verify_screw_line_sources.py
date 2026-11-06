#!/usr/bin/env python3
"""Exact bounded source certificate for the eight P21 scalar beta physical pairs.

Promoted from the independent R61 native contract probe. This is raw archive
evidence, not a production subduction or all-embedding certificate. Run under
prlimit --as=536870912:536870912; --inventory emits decoded native payloads.
The full ZIP hashes are pinned, but only the strict SG1..4 prefix is materialized.
"""
import hashlib, io, json, resource, sys, zipfile
from collections import Counter
from fractions import Fraction as F
from itertools import product
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import generate_irrep_data as g
from scripts.iso_irrep_exact import SourceArchive, parse_exact_source_lines
from scripts.spinor_exact import Complex24 as C, Radical24 as R, ZETA24_POWERS
assert resource.getrlimit(resource.RLIMIT_AS) == (536870912, 536870912)
ZERO, ONE, I = C(), C(R(F(1))), C(im=R(F(1)))
I2 = ((ONE,ZERO),(ZERO,ONE))
counts=Counter()
def phase(q):
    n=(q%1)*24
    assert n.denominator==1, q
    return ZETA24_POWERS[int(n)]
def radical(v): return R(*(F(x,4) for x in v))
def mm(a,b): return tuple(tuple(sum((a[i][k]*b[k][j] for k in range(len(b))),ZERO) for j in range(len(b[0]))) for i in range(len(a)))
def adj(a): return tuple(tuple(v.conjugate() for v in col) for col in zip(*a))
def add(a,b): return tuple(x+y for x,y in zip(a,b))
def sub(a,b): return tuple(x-y for x,y in zip(a,b))
def mv(a,b): return tuple(sum(x*y for x,y in zip(row,b)) for row in a)
def exact_matrices(r, lines):
    idx=next(i for i,l in enumerate(lines) if '"' in l and int(l.split()[0])==r.irnumber)
    if r.archive is SourceArchive.PIR:
        raw,idx=g._read_exact_pir_int_block(lines,idx+1,16*len(r.k_arms),'native k')
    else:
        raw,idx=g._read_exact_cir_block(lines,idx+1,16*len(r.k_arms),'native k',g._parse_cir_integer)
    assert tuple(raw)==tuple(x for a in r.k_arms for x in a.raw_augmented)
    result=[]
    for op,irt in zip(r.operations,r.irtranslations):
        if r.archive is SourceArchive.PIR:
            ro,rt,_,sp,idx=g._read_pir_operation_payload(lines,idx,r.dimension,r.special,'native matrix')
            vals=tuple(C(radical(g._decode_pir_matrix_token(s))) for s in sp)
        else:
            ro,idx=g._parse_cir_operation_row(lines,idx,'native operation')
            rt=None
            if not r.special: rt,idx=g._parse_cir_irtranslation_row(lines,idx,'native irtranslation')
            raw,_,idx=g._parse_cir_complex_block(lines,idx,r.dimension**2,'native matrix')
            vals=tuple(C(radical(re),radical(im)) for re,im in raw)
        assert tuple(ro)==op.raw_augmented
        assert (None if rt is None else tuple(rt))==(None if irt is None else irt.raw)
        result.append(tuple(tuple(vals[i*r.dimension+j] for j in range(r.dimension)) for i in range(r.dimension)))
        counts['native_operation_bindings']+=1
        counts['native_matrix_tokens']+=len(vals)
    return tuple(result)
def prefix(archive):
    name=archive.value+'_data.zip'; data=(ROOT/'isotropy_subgroup'/name).read_bytes()
    assert hashlib.sha256(data).hexdigest()==g.PINNED_ARCHIVE_SHA256[name]
    lines=[]
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        assert z.namelist().count(archive.value+'_data.txt')==1
        with z.open(archive.value+'_data.txt') as f:
            for raw in f:
                assert raw.endswith(b'\n')
                line=raw.decode('ascii').rstrip('\r\n')
                if '"' in line:
                    sg=(g._parse_pir_header(line,len(lines)+1)[1] if archive is SourceArchive.PIR else g._parse_cir_header(line,len(lines)+1)['sg'])
                    if sg>4: break
                lines.append(line)
                assert len(lines)<5000
    return lines,parse_exact_source_lines(lines,archive,validate_census=False)
DB={}; LINES={}; MATRICES={}
for archive in SourceArchive:
    lines,records=prefix(archive); LINES[archive]=lines; DB[archive]=records
    assert records[0].irnumber==1 and max(r.spacegroup for r in records)==4
    for r in records:
        if r.spacegroup==4:
            MATRICES[archive,r.irnumber]=exact_matrices(r,lines)
# Explicit native identities: no inference from PIR=CIR numbering.
SPECS=(('LD1LE1','LD1',56,56,0,0,1),('LD2LE2','LD2',57,57,0,0,-1),
       ('U1UA1','U1',68,72,1,1,1),('U2UA2','U2',69,73,1,1,-1),
       ('V1VA1','V1',70,74,0,1,1),('V2VA2','V2',71,75,0,1,-1),
       ('W1WA1','W1',72,76,1,0,1),('W2WA2','W2',73,77,1,0,-1))
pirs={r.irnumber:r for r in DB[SourceArchive.PIR] if r.spacegroup==4}
cirs={r.irnumber:r for r in DB[SourceArchive.CIR] if r.spacegroup==4}
ops=pirs[56].operations
assert tuple(o.rotation for o in ops)==(((1,0,0),(0,1,0),(0,0,1)),((-1,0,0),(0,1,0),(0,0,-1)))
assert tuple(o.translation for o in ops)==((F(0),)*3,(F(0),F(1,2),F(0)))
for a,records,ids in [(SourceArchive.PIR,pirs,{s[2] for s in SPECS}),(SourceArchive.CIR,cirs,{s[3] for s in SPECS})]:
    actual={r.irnumber for r in records.values() if not r.special and r.irtype==3 and r.dimension==(2 if a is SourceArchive.PIR else 1) and len(r.k_arms)==1 and r.k_arms[0].parameters==(None,(F(0),F(1),F(0)),None)}
    assert actual==ids
for pn,cn,pi,ci,x,z,sign in SPECS:
    for r,name,dim in [(pirs[pi],pn,2),(cirs[ci],cn,1)]:
        assert r.irrep_label==name and r.dimension==dim and r.centering.value=='P'
        assert r.irtype==3 and r.kcount==r.pmkcount==len(r.k_arms)==1
        assert r.operations==ops and all(t.vector==o.translation for t,o in zip(r.irtranslations,ops))
        assert r.k_arms[0].constant==(F(x,2),F(0),F(z,2))
        identity=I2 if dim==2 else ((ONE,),)
        assert MATRICES[r.archive,r.irnumber]==(identity,tuple(tuple(sign*v for v in row) for row in identity))
    assert pn[len(cn):] not in {r.irrep_label for r in cirs.values()}
def kval(r,b):
    a=r.k_arms[0]
    return tuple(a.constant[i]+(a.parameters[1][i]*b if a.parameters[1] else 0) for i in range(3))
def rep(r,b,op):
    rot,t=op; j=next(j for j,o in enumerate(r.operations) if o.rotation==rot)
    delta=sub(t,r.operations[j].translation)
    if any(F(x).denominator!=1 for x in delta): raise ValueError('non-P-lattice operation')
    if r.irtranslations[j] is not None: delta=add(delta,r.irtranslations[j].vector)
    z=phase(sum(x*y for x,y in zip(kval(r,b),delta)))
    trans=((C(z.re),C(z.im)),(C(-z.im),C(z.re))) if r.archive is SourceArchive.PIR and r.dimension==2 else ((z,),)
    return mm(trans,MATRICES[r.archive,r.irnumber][j])
def compose(a,b):
    ra,ta=a; rb,tb=b
    rot=tuple(tuple(sum(ra[i][k]*rb[k][j] for k in range(3)) for j in range(3)) for i in range(3))
    return rot,add(ta,mv(ra,tb))
s=C(R(b=F(1,2))); U=((s,-I*s),(s,I*s)); Q=((ONE,ZERO),(ZERO,-ONE))
assert mm(U,adj(U))==I2
translations=((0,0,0),(1,0,0),(0,1,0),(0,0,1),(-1,2,-1))
operations=tuple((o.rotation,add(o.translation,t)) for o,t in product(ops,translations))
parameters=tuple(map(F,['0','1/6','1/4','1/3','1/2','3/4','1','3/2','-1/4']))
endpoints=[]; norms=[]
if '--inventory' not in sys.argv:
    for pn,cn,pi,ci,xbit,zbit,sign in SPECS:
        p,c=pirs[pi],cirs[ci]; peer=cirs[ci+sign]
        for b in parameters:
            # Per-model finite cache, discarded after each beta; never full-db/unbounded.
            cache={}
            def evaluate(r,op):
                key=(r.archive,r.irnumber,op)
                if key not in cache: cache[key]=rep(r,b,op)
                assert len(cache)<1000
                return cache[key]
            k=kval(c,b)
            assert len({tuple(v%1 for v in mv(o.rotation,k)) for o in ops})==1
            assert sum(all(v.denominator==1 for v in sub(mv(o.rotation,k),k)) for o in ops)==2
            assert len({tuple(v%1 for v in k),tuple((-v)%1 for v in k)})==(1 if (2*b).denominator==1 else 2)
            counts['little_group_and_arm_models']+=1
            for op in operations:
                P=evaluate(p,op); z=evaluate(c,op)[0][0]
                assert mm(P,adj(P))==I2
                assert mm(U,mm(P,adj(U)))==((z,ZERO),(ZERO,z.conjugate()))
                assert P[0][0]+P[1][1]==z+z.conjugate()
                assert rep(c,-b,op)[0][0]==z.conjugate()
                assert mm(Q,mm(P,Q))==rep(p,-b,op)
                counts['complexification_trace_conjugacy']+=1
                assert rep(c,b+2,op)==evaluate(c,op)
                assert rep(c,b+1,op)==rep(peer,b,op)
                counts['period_two_and_label_swap']+=1
            for a,d in product(operations,repeat=2):
                ad=compose(a,d)
                assert mm(evaluate(p,a),evaluate(p,d))==evaluate(p,ad)
                za,zd,zad=(evaluate(c,o)[0][0] for o in (a,d,ad))
                assert za*zd==zad and za.conjugate()*zd.conjugate()==zad.conjugate()
                counts['pir_group_products']+=1; counts['cir_conjugate_group_products']+=2
            finite=tuple((o.rotation,add(o.translation,(F(x),F(y),F(z)))) for o,x,y,z in product(ops,range(1+xbit),range(b.denominator),range(1+zbit)))
            normp=normc=ZERO
            for op in finite:
                z=rep(c,b,op)[0][0]; chi=z+z.conjugate()
                normp+=chi*chi.conjugate(); normc+=z*z.conjugate()
            expected=4 if b.denominator==1 else 2
            assert normp==len(finite)*expected*ONE and normc==len(finite)*ONE
            counts['finite_quotient_norm_models']+=1
            norms.append([pi,str(b),len(finite),expected])
            for bad in [(ops[0].rotation,(F(0),F(1,2),F(0))),(ops[1].rotation,(F(0),)*3)]:
                try: rep(c,b,bad)
                except ValueError: counts['nonlattice_rejections']+=1
                else: raise AssertionError('nonlattice accepted')
        for b in map(F,['0','1/2','1','3/2','2']):
            k=kval(c,b); rows=[]
            for conjugate in (False,True):
                hits=[]
                for e in cirs.values():
                    if not e.special or e.dimension!=1 or any(v.denominator!=1 for v in sub(k,e.k_arms[0].constant)): continue
                    if all(rep(e,F(0),op)[0][0]==(rep(c,b,op)[0][0].conjugate() if conjugate else rep(c,b,op)[0][0]) for op in operations): hits.append(e)
                assert len(hits)==1
                e=hits[0]
                owners=[o for o in pirs.values() if o.special and o.k_arms==e.k_arms and (o.irrep_label==e.irrep_label or o.irrep_label in ('C1C2','D1D2','E1E2','Z1Z2') and e.irrep_label in o.irrep_label)]
                assert len(owners)==1
                owner=owners[0]
                rows.append([e.irnumber,e.irrep_label,owner.irnumber,owner.irrep_label,None if owner.dimension==1 else (0 if owner.irrep_label.startswith(e.irrep_label) else 1)])
                for op in operations:
                    O=rep(owner,F(0),op)
                    if owner.dimension==1:
                        assert O==rep(e,F(0),op)
                    else:
                        oz=rep(c,b,op)[0][0]
                        OU=mm(U,mm(O,adj(U)))
                        assert OU in (((oz,ZERO),(ZERO,oz.conjugate())),((oz.conjugate(),ZERO),(ZERO,oz)))
                    counts['native_owner_matrix_checks']+=1
                counts['parent_endpoint_component_bindings']+=1
            assert (rows[0]==rows[1])==(b.denominator==1)
            endpoints.append({'pir':pi,'cir':ci,'b':str(b),'seed':rows[0],'conjugate':rows[1]})
        screw=(ops[1].rotation,ops[1].translation); ty=compose(screw,screw)
        assert rep(c,F(1,2),screw)[0][0]==sign*I
        assert rep(c,F(1,2),ty)[0][0]==-ONE
        assert (sign*ONE)*(sign*ONE)!=rep(c,F(1,2),ty)[0][0]
        counts['omit_irtranslation_group_law_controls']+=1
        assert rep(c,F(0),screw)!=rep(c,F(1),screw)
        assert rep(p,F(0),screw)!=rep(p,F(1),screw)
        counts['false_period_one_controls']+=1
        assert 2*rep(c,F(1,2),screw)[0][0]!=sum((rep(p,F(1,2),screw)[i][i] for i in range(2)),ZERO)
        counts['false_double_CIR_controls']+=1
# Exact representability caution and overflow-free split comparison witnesses.
MAX=2**127-1; MIN=-2**127
assert (F(-1,MAX)%2).numerator>MAX
extremes=[F(n,d) for n in [MIN,MIN+1,-MAX,-1,0,1,MAX-1,MAX] for d in [1,2,3,2**90,2**120,MAX-1,MAX]]
near=[F(i,4)+F(s,2**power) for i in range(-8,9) for s in [-1,0,1] for power in [90,120]]
for b in extremes+near:
    q,r=divmod(b.numerator,b.denominator)
    key=(q%2,F(r,b.denominator))
    exact=b%2; assert key==(exact.numerator//exact.denominator,exact%1)
    for cut in map(F,['0','1/4','1/2','3/4','1','5/4','3/2','7/4']):
        cq,cr=divmod(cut.numerator,cut.denominator)
        assert (key==(cq,F(cr,cut.denominator)))==(exact==cut)
        counts['split_mod2_exact_comparisons']+=1

def scalar_json(c): return [list(map(str,(c.re.a,c.re.b,c.re.c,c.re.d))),list(map(str,(c.im.a,c.im.b,c.im.c,c.im.d)))]
def matrix_json(m): return [[scalar_json(v) for v in row] for row in m]
report={'scope':'pinned SG1..4 strict prefix; all SG4 native point payloads; no production execution','prefix_record_counts':{a.value:len(rs) for a,rs in DB.items()},'prefix_line_counts':{a.value:len(ls) for a,ls in LINES.items()},'records':[{'archive':r.archive.value,'id':r.irnumber,'label':r.irrep_label,'dimension':r.dimension,'special':r.special,'irtype':r.irtype,'kcount':r.kcount,'pmkcount':r.pmkcount,'k0':list(map(str,r.k_arms[0].constant)),'parameters':[[str(x) for x in v] if v is not None else None for v in r.k_arms[0].parameters],'operations':[{'rotation':o.rotation,'translation':list(map(str,o.translation)),'irtranslation':list(map(str,t.vector)) if t else None,'matrix':matrix_json(m)} for o,t,m in zip(r.operations,r.irtranslations,MATRICES[r.archive,r.irnumber])]} for rs in DB.values() for r in rs if r.spacegroup==4],'counts':dict(counts),'maxrss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
if '--inventory' not in sys.argv:
    expected = {'pir_group_products':7200, 'cir_conjugate_group_products':14400,
        'complexification_trace_conjugacy':720, 'period_two_and_label_swap':720,
        'finite_quotient_norm_models':72, 'parent_endpoint_component_bindings':80,
        'native_owner_matrix_checks':800, 'nonlattice_rejections':144,
        'omit_irtranslation_group_law_controls':8, 'false_period_one_controls':8,
        'false_double_CIR_controls':8, 'split_mod2_exact_comparisons':1264}
    for key, value in expected.items():
        assert counts[key] == value, (key, counts[key], value)
    report.pop('records')
    report.update({'status':'PASS','source_specs':SPECS,'sample_parameters':list(map(str,parameters)),'operations_per_model':len(operations),'counts':dict(counts),'endpoint_bindings':endpoints,'norm_models':norms,'arithmetic_unrepresentable_example':str(F(-1,MAX)%2),'input_sha256':{str(p):hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['scripts/iso_irrep_exact.py','scripts/generate_irrep_data.py','scripts/spinor_exact.py','isotropy_subgroup/PIR_data.zip','isotropy_subgroup/CIR_data.zip']},'maxrss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
print(json.dumps(report,indent=2))
