#!/usr/bin/env python3
"""Exact native proof for all eight SG3 scalar physical conjugate lines.

Uses pinned archive readers and existing Q(zeta_24) arithmetic; no runtime iso
call or new solver. Pair with tests/physical_line_subduction.rs to bind these
authoritative source fields to the immutable Rust registry. This script alone
is source evidence, not target-subduction or complete parameter-domain coverage.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.iso_irrep_exact import load_exact_iso_irrep_sources
from scripts.audit_irrep_archives import rank
from scripts import generate_irrep_data as g
from scripts.spinor_exact import Complex24 as C, Radical24 as R, ZETA24_POWERS

ZERO, ONE, I = C(), C(R(F(1))), C(im=R(F(1)))
I2 = ((ONE, ZERO), (ZERO, ONE))
counts = Counter()

def addv(a, b): return tuple(x+y for x,y in zip(a,b))
def subv(a, b): return tuple(x-y for x,y in zip(a,b))
def mv(a, v): return tuple(sum(x*y for x,y in zip(row,v)) for row in a)
def mm(a, b):
    return tuple(tuple(sum((a[i][k]*b[k][j] for k in range(len(b))), ZERO)
                       for j in range(len(b[0]))) for i in range(len(a)))
def adj(a): return tuple(tuple(v.conjugate() for v in col) for col in zip(*a))
def tr(a): return sum((a[i][i] for i in range(len(a))), ZERO)
def phase(q):
    n = (q % 1)*24
    assert n.denominator == 1, q
    return ZETA24_POWERS[int(n)]
def radical(v): return R(*(F(x,4) for x in v))
def exact_record_matrices(record, lines):
    # Source model intentionally omits matrices: reuse pinned generator readers.
    idx = next(i for i,l in enumerate(lines) if '"' in l and int(l.split()[0]) == record.irnumber)
    if record.archive.value == 'PIR':
        raw, idx = g._read_exact_pir_int_block(lines, idx+1, 16*len(record.k_arms), 'native k')
    else:
        raw, idx = g._read_exact_cir_block(lines, idx+1, 16*len(record.k_arms), 'native k', g._parse_cir_integer)
    assert tuple(raw) == tuple(x for a in record.k_arms for x in a.raw_augmented)
    matrices = []
    for op, irt in zip(record.operations, record.irtranslations):
        if record.archive.value == 'PIR':
            rawop, rawt, _, spellings, idx = g._read_pir_operation_payload(lines, idx, record.dimension, record.special, 'native matrix')
            vals = tuple(C(radical(g._decode_pir_matrix_token(s))) for s in spellings)
        else:
            rawop, idx = g._parse_cir_operation_row(lines, idx, 'native operation')
            rawt = None
            if not record.special: rawt, idx = g._parse_cir_irtranslation_row(lines, idx, 'native irtranslation')
            rawvals, _, idx = g._parse_cir_complex_block(lines, idx, record.dimension**2, 'native matrix')
            vals = tuple(C(radical(re),radical(im)) for re,im in rawvals)
        assert tuple(rawop) == op.raw_augmented
        assert (None if rawt is None else tuple(rawt)) == (None if irt is None else irt.raw)
        matrices.append(tuple(tuple(vals[i*record.dimension+j] for j in range(record.dimension)) for i in range(record.dimension)))
        counts['raw_operation_bindings'] += 1
        counts['exact_matrix_tokens'] += len(vals)
    return tuple(matrices)

def evaluate_k(r, b):
    a = r.k_arms[0]
    return tuple(a.constant[i]+sum((v[i]*p for v,p in zip(a.parameters,(F(0),b,F(0))) if v is not None), F(0)) for i in range(3))

translations = ((0,0,0),(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1),(1,1,1),(1,2,-1))
parameters = tuple(map(F, ['1/3','-1/3','1/4','-1/4','0','1/2','1']))
db = load_exact_iso_irrep_sources()
pirs = {r.irrep_label:r for r in db.pir_records if r.spacegroup == 3}
cirs = {r.irrep_label:r for r in db.cir_records if r.spacegroup == 3}
# Explicit native identities, transverse half-offset bits, and C2 point sign.
SOURCES = (
    ('LD1LE1','LD1',29,0,0,1), ('LD2LE2','LD2',30,0,0,-1),
    ('U1UA1','U1',45,1,1,1), ('U2UA2','U2',46,1,1,-1),
    ('V1VA1','V1',47,0,1,1), ('V2VA2','V2',48,0,1,-1),
    ('W1WA1','W1',49,1,0,1), ('W2WA2','W2',50,1,0,-1),
)
specs = {irn:(x,z,sign) for _,_,irn,x,z,sign in SOURCES}
# Exhaust this *scalar, single-arm beta-line* cohort, not all SG3 parameters.
for archive, records, names in [('PIR',pirs.values(),{s[0] for s in SOURCES}),('CIR',cirs.values(),{s[1] for s in SOURCES})]:
    actual = {r.irrep_label for r in records if not r.special and r.irtype==3
              and r.dimension==(2 if archive=='PIR' else 1) and len(r.k_arms)==1
              and r.k_arms[0].parameters==(None,(F(0),F(1),F(0)),None)}
    assert actual==names
for pname,cname,irn,x,z,sign in SOURCES:
    assert (pirs[pname].irnumber,cirs[cname].irnumber)==(irn,irn)
    assert pirs[pname].k_arms==cirs[cname].k_arms
    assert pname[len(cname):] not in cirs
selected = [pirs[s[0]] for s in SOURCES]+[cirs[s[1]] for s in SOURCES]
selected += [r for r in [*pirs.values(),*cirs.values()] if r.special]
pir_lines, cir_lines = g._read_pir_lines(), g._read_cir_lines()
matrices = {(r.archive.value,r.irnumber):exact_record_matrices(r,pir_lines if r.archive.value=='PIR' else cir_lines) for r in selected}
ops = pirs['LD1LE1'].operations
assert len(ops)==2
assert ops[0].rotation==((1,0,0),(0,1,0),(0,0,1))
assert ops[1].rotation==((-1,0,0),(0,1,0),(0,0,-1))
assert all(o.translation==(F(0),)*3 for o in ops)
for r in selected:
    assert r.operations==ops and len(r.k_arms)==r.kcount==r.pmkcount==1
    expected = I2 if r.dimension==2 else ((ONE,),)
    sign = specs[r.irnumber][2] if not r.special else (-1 if r.irrep_label.endswith("2") else 1)
    assert matrices[r.archive.value,r.irnumber]==(expected,tuple(tuple(sign*v for v in row) for row in expected))
    if not r.special:
        assert r.centering.value == 'P'
        assert r.dimension == (2 if r.archive.value == 'PIR' else 1)
        x,z,_ = specs[r.irnumber]
        assert r.k_arms[0].constant == (F(x,2),F(0),F(z,2))
        assert r.k_arms[0].parameters == (None,(F(0),F(1),F(0)),None)
        assert r.irtype==3 and rank(r.k_arms[0].parameters)==1
        assert all(t.vector==(F(0),)*3 for t in r.irtranslations)

@lru_cache(None)
def rep(archive, irn, b, rot, t):
    r = (db.pir_records if archive=='PIR' else db.cir_records)[irn-1]
    i = next(i for i,o in enumerate(r.operations) if o.rotation==rot)
    delta = subv(t,r.operations[i].translation)
    if any(x.denominator!=1 for x in map(F,delta)):
        raise ValueError('not a P-lattice operation')
    if r.irtranslations[i] is not None: delta=addv(delta,r.irtranslations[i].vector)
    z = phase(sum(x*y for x,y in zip(evaluate_k(r,b),delta)))
    trans = ((C(z.re),C(z.im)),(C(-z.im),C(z.re))) if archive=='PIR' and r.dimension==2 else ((z,),)
    return mm(trans, matrices[archive,irn][i])

operations = tuple((o.rotation,tuple(map(F,t))) for o,t in product(ops,translations))
def compose(a,b):
    ra,ta=a; rb,tb=b
    rot=tuple(tuple(sum(ra[i][k]*rb[k][j] for k in range(3)) for j in range(3)) for i in range(3))
    return rot,addv(ta,mv(ra,tb))
# Official PIR->complex map U, first eigenvalue +phase, second -phase.
s = C(R(b=F(1,2)))
U=((s,-I*s),(s,I*s))
assert mm(U,adj(U))==I2
samples=[]
for pname,cname,irnumber,xbit,zbit,point_sign in SOURCES:
    p,c=pirs[pname],cirs[cname]
    for b in parameters:
        pk=evaluate_k(p,b)
        images={mv(o.rotation,pk) for o in ops} # here R^-T == R
        modimages={tuple(x%1 for x in a) for a in images}
        physical={tuple(x%1 for x in pk),tuple((-x)%1 for x in pk)}
        assert len(modimages)==1
        assert len(physical)==(1 if (2*b).denominator==1 else 2)
        little=sum(all(x.denominator==1 for x in subv(mv(o.rotation,pk),pk)) for o in ops)
        assert little==2
        counts['little_group_and_physical_arm_samples']+=1
        for rot,t in operations:
            P=rep('PIR',p.irnumber,b,rot,t); z=rep('CIR',c.irnumber,b,rot,t)[0][0]
            assert mm(P,adj(P))==I2
            assert mm(U,mm(P,adj(U)))==((z,ZERO),(ZERO,z.conjugate()))
            counts['exact_complexification_matches']+=1
            assert tr(P)==z+z.conjugate()
            counts['exact_character_pair_matches']+=1
            assert rep('CIR',c.irnumber,-b,rot,t)[0][0]==z.conjugate()
            counts['negative_parameter_conjugacy']+=1
            if (2*b).denominator==1:
                endpoints=[e for e in cirs.values() if e.special
                           and e.k_arms[0].constant==tuple(x%1 for x in pk)
                           and matrices['CIR',e.irnumber][1]==((point_sign*ONE,),)]
                assert len(endpoints)==1
                e=endpoints[0]
                w=rep('CIR',e.irnumber,F(0),rot,t)[0][0]
                assert z==w and P==((w,ZERO),(ZERO,w))
                counts['endpoint_duplicate_real_irrep_matches']+=1
        for a,d in product(operations,repeat=2):
            ad=compose(a,d)
            assert mm(rep('PIR',p.irnumber,b,*a),rep('PIR',p.irnumber,b,*d))==rep('PIR',p.irnumber,b,*ad)
            counts['exact_pir_group_products']+=1
            za=rep('CIR',c.irnumber,b,*a)[0][0]; zd=rep('CIR',c.irnumber,b,*d)[0][0]; zad=rep('CIR',c.irnumber,b,*ad)[0][0]
            assert za*zd==zad
            assert za.conjugate()*zd.conjugate()==zad.conjugate()
            counts['exact_cir_or_conjugate_group_products']+=2
        norm_p,norm_c=ZERO,ZERO
        finite_ops=tuple((o.rotation,(F(x),F(y),F(z))) for o,x,y,z in product(ops,range(1+xbit),range(b.denominator),range(1+zbit)))
        for rot,t in finite_ops:
            chi=tr(rep('PIR',p.irnumber,b,rot,t)); z=rep('CIR',c.irnumber,b,rot,t)[0][0]
            norm_p+=chi*chi.conjugate(); norm_c+=z*z.conjugate()
            counts['exact_character_norm_terms']+=1
        expected_norm=4 if (2*b).denominator==1 else 2
        assert norm_p==len(finite_ops)*expected_norm*ONE and norm_c==len(finite_ops)*ONE
        counts['exact_character_norm_models']+=1
        try: rep('PIR',p.irnumber,b,ops[0].rotation,(F(0),F(1,2),F(0)))
        except ValueError: counts['nonlattice_translation_rejections']+=1
        else: raise AssertionError('nonlattice operation accepted')
        # Counterfactual dim2 same-arm model agrees on E,C2 but fails generic Ty.
        z=rep('CIR',c.irnumber,b,ops[0].rotation,(F(0),F(1),F(0)))[0][0]
        realtrace=tr(rep('PIR',p.irnumber,b,ops[0].rotation,(F(0),F(1),F(0))))
        assert (2*z==realtrace)==((2*b).denominator==1)
        counts['false_doubled_cir_detected' if 2*z!=realtrace else 'special_parameters_hide_doubled_cir_error']+=1
        if (2*b).denominator!=1:
            # Synthetic re-section C2: stored tau=(0,1,0), irtranslation=tau.
            # Correct t-tau+irtranslation=t; dropping irtranslation is wrong.
            assert phase(b*(-1+1))==ONE and phase(-b)!=ONE
            counts['synthetic_nonzero_irtranslation_omission_detected']+=1
        # Equivalent +/-b real matrices via fixed reflection Q=diag(1,-1).
        Q=((ONE,ZERO),(ZERO,-ONE))
        assert mm(Q,mm(rep('PIR',p.irnumber,b,ops[0].rotation,(F(0),F(1),F(0))),Q))==rep('PIR',p.irnumber,-b,ops[0].rotation,(F(0),F(1),F(0)))
        samples.append({'pir':p.irnumber,'cir':c.irnumber,'b':str(b),'physical_points':len(physical),'parent_orbit_points':len(modimages),'little_order':little,'complex_character_norm_pir':expected_norm,'complex_character_norm_cir':1,'Ty_trace':str(realtrace)})
# Distinguish the actual two independent offset bits AND the C2 point sign.
for _,_,irn,x,z,sign in SOURCES:
    for b in parameters:
        for t,bit in [((F(1),F(0),F(0)),x),((F(0),F(0),F(1)),z)]:
            assert rep('CIR',irn,b,ops[0].rotation,t)[0][0]==(1-2*bit)*ONE
            counts['offset_axis_sign_checks']+=1
        assert rep('CIR',irn,b,ops[1].rotation,(F(0),)*3)[0][0]==sign*ONE
        counts['point_sign_checks']+=1
        if sign==-1:
            assert rep('CIR',irn,b,ops[1].rotation,(F(0),)*3)!=rep('CIR',irn-1,b,ops[1].rotation,(F(0),)*3)
            counts['wrong_positive_point_source_detected']+=1
report={'scope':'all eight SG3 scalar beta physical conjugate lines; not general affine or all continuous SG3 sources','source_specs':[{'pir':p,'cir':c,'id':i,'offset_bits':[x,z],'C2_sign':s} for p,c,i,x,z,s in SOURCES], 'source_census':[len(db.pir_records),len(db.cir_records)],'sg3_cir_labels':[[r.irnumber,r.irrep_label] for r in cirs.values()], 'selected_records':[{'archive':r.archive.value,'id':r.irnumber,'label':r.irrep_label,'dim':r.dimension,'irtype':r.irtype,'kcount':r.kcount,'pmkcount':r.pmkcount,'arms':len(r.k_arms),'constant':list(map(str,r.k_arms[0].constant)),'parameters':[[str(v) for v in x] if x is not None else None for x in r.k_arms[0].parameters],'irtranslations':[list(map(str,t.vector)) if t else None for t in r.irtranslations]} for r in selected], 'sample_parameters':list(map(str,parameters)),'translation_count':len(translations),'operation_count_per_model':len(operations),'counts':dict(counts),'samples':samples,'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ('scripts/iso_irrep_exact.py','scripts/audit_irrep_archives.py','scripts/generate_irrep_data.py','scripts/spinor_exact.py','isotropy_subgroup/PIR_data.zip','isotropy_subgroup/CIR_data.zip')},'limitations':'Native SG3 sources only; finite exact witnesses, not production subduction or full battery. Existing Complex24 used only as exact arithmetic, not spinor evidence.'}
print(json.dumps(report,indent=2))
