#!/usr/bin/env python3
"""Exact bounded native certificate for four C-centered scalar beta pairs.

Pins whole ZIP identities and parses only the strict SG1..5 prefix. Requires
AS512MiB. The native point matrices and complete centered-translation phase
are checked independently of Rust; this is not a child-endpoint or production
subduction certificate, nor a full-archive census rerun.
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
                    if sg>5: break
                lines.append(line)
                assert len(lines)<5000
    return lines,parse_exact_source_lines(lines,archive,validate_census=False)

SPECS=(('LD1LE1','LD1',79,83,0,1),('LD2LE2','LD2',80,84,0,-1),('U1UA1','U1',89,93,1,1),('U2UA2','U2',90,94,1,-1))
DB={};MATRICES={}
for archive in SourceArchive:
    lines,records=prefix(archive);DB[archive]={r.irnumber:r for r in records if r.spacegroup==5}
    ids={s[2] if archive is SourceArchive.PIR else s[3] for s in SPECS}
    actual={r.irnumber for r in DB[archive].values() if not r.special and r.irtype==3 and r.dimension==(2 if archive is SourceArchive.PIR else 1) and len(r.k_arms)==1 and r.k_arms[0].parameters==(None,(F(0),F(1),F(0)),None)}
    assert actual==ids
    for irn in ids:MATRICES[archive,irn]=exact_matrices(DB[archive][irn],lines)
pirs,cirs=DB[SourceArchive.PIR],DB[SourceArchive.CIR];ops=pirs[79].operations
assert tuple(o.rotation for o in ops)==(((1,0,0),(0,1,0),(0,0,1)),((-1,0,0),(0,1,0),(0,0,-1)))
assert all(o.translation==(F(0),)*3 for o in ops)
for pn,cn,pi,ci,zbit,sign in SPECS:
    for r,name,dim in [(pirs[pi],pn,2),(cirs[ci],cn,1)]:
        assert r.irrep_label==name and r.dimension==dim and r.centering.value=='C'
        assert r.irtype==3 and r.kcount==r.pmkcount==len(r.k_arms)==1
        assert r.operations==ops and all(t.vector==(F(0),)*3 for t in r.irtranslations)
        assert r.k_arms[0].constant==(F(0),F(0),F(zbit,2))
        identity=I2 if dim==2 else ((ONE,),)
        assert MATRICES[r.archive,r.irnumber]==(identity,tuple(tuple(sign*v for v in row) for row in identity))
    assert pn[len(cn):] not in {r.irrep_label for r in cirs.values()}
def lattice(t):
    x,y,z=map(F,t)
    return (2*x).denominator==(2*y).denominator==z.denominator==1 and (x-y).denominator==1
def reciprocal(k):return all(F(v).denominator==1 for v in k) and (k[0]+k[1])%2==0
def kval(r,b):return add(r.k_arms[0].constant,(F(0),b,F(0)))
def rep(r,b,op):
    rot,t=op;j=next(i for i,o in enumerate(r.operations) if o.rotation==rot)
    delta=sub(t,r.operations[j].translation)
    if not lattice(delta):raise ValueError('non-C-lattice operation')
    delta=add(delta,r.irtranslations[j].vector);z=phase(sum(x*y for x,y in zip(kval(r,b),delta)))
    trans=((C(z.re),C(z.im)),(C(-z.im),C(z.re))) if r.archive is SourceArchive.PIR else ((z,),)
    return mm(trans,MATRICES[r.archive,r.irnumber][j])
def compose(a,b):
    ra,ta=a;rb,tb=b
    return tuple(tuple(sum(ra[i][k]*rb[k][j] for k in range(3)) for j in range(3)) for i in range(3)),add(ta,mv(ra,tb))
s=C(R(b=F(1,2)));U=((s,-I*s),(s,I*s));assert mm(U,adj(U))==I2
center=(F(1,2),F(1,2),F(0));assert lattice(center) and not reciprocal((F(0),F(1),F(0))) and reciprocal((F(0),F(2),F(0)))
translations=((0,0,0),(1,0,0),(0,1,0),(0,0,1),center,(F(-3,2),F(5,2),F(-1)))
operations=tuple((o.rotation,tuple(map(F,t))) for o,t in product(ops,translations))
parameters=tuple(map(F,['0','1/6','1/4','1/3','1/2','3/4','1','3/2','-1/4']))
for pn,cn,pi,ci,zbit,sign in SPECS:
    p,c=pirs[pi],cirs[ci];peer=cirs[ci+sign]
    for b in parameters:
        cache={}
        def evaluate(r,op):
            key=(r.archive,r.irnumber,op)
            if key not in cache:cache[key]=rep(r,b,op)
            assert len(cache)<1000
            return cache[key]
        k=kval(c,b);assert all(reciprocal(sub(mv(o.rotation,k),k)) for o in ops)
        assert reciprocal(tuple(2*v for v in k))==(b.denominator==1)
        counts['parent_little_group_and_collision_models']+=1
        for op in operations:
            P=evaluate(p,op);z=evaluate(c,op)[0][0]
            assert mm(P,adj(P))==I2 and mm(U,mm(P,adj(U)))==((z,ZERO),(ZERO,z.conjugate()))
            assert P[0][0]+P[1][1]==z+z.conjugate()
            assert rep(c,-b,op)[0][0]==z.conjugate() and rep(c,b+2,op)[0][0]==z
            counts['complexification_and_period2_checks']+=1
        for a,other in product(operations,repeat=2):
            ab=compose(a,other);assert lattice(ab[1])
            for r in [p,c]:
                assert mm(evaluate(r,a),evaluate(r,other))==evaluate(r,ab)
                counts['PIR_products' if r is p else 'CIR_products']+=1
        center_op=(ops[0].rotation,center);z=rep(c,b,center_op)[0][0]
        assert rep(c,b+1,center_op)[0][0]==-z!=z
        assert rep(c,b+1,center_op)[0][0]!=rep(peer,b,center_op)[0][0]
        counts['period1_and_sign_swap_negative_controls']+=1
        for bad in [(F(1,2),0,0),(0,F(1,2),0),(0,0,F(1,2))]:
            for r in [p,c]:
                try:rep(r,b,(ops[0].rotation,bad))
                except ValueError:counts['nonlattice_rejections']+=1
                else:raise AssertionError('accepted non-C translation')
        assert rep(c,b,(ops[1].rotation,(F(0),)*3))[0][0]==sign*ONE
        counts['point_sign_checks']+=1
        if zbit:
            assert rep(c,b,(ops[0].rotation,(F(0),F(0),F(1))))[0][0]==-ONE
            counts['omitted_offset_negative_controls']+=1
        # Finite quotient of the source kernel: integer x is trivial; C2 acts trivially on its quotient.
        N=(b/2).denominator;norm_c=ZERO;norm_p=ZERO
        assert rep(c,b,(ops[0].rotation,(F(1),F(0),F(0))))[0][0]==ONE
        for o,n,q in product(ops,range(N),range(2)):
            z=rep(c,b,(o.rotation,(F(n,2),F(n,2),F(q))))[0][0];chi=z+z.conjugate()
            norm_c+=z*z.conjugate();norm_p+=chi*chi.conjugate()
        assert norm_c==C(R(F(4*N)))
        assert norm_p==C(R(F(4*N*(4 if b.denominator==1 else 2))))
        counts['finite_quotient_norms']+=1
pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__).resolve(),ROOT/'scripts/verify_screw_line_sources.py',ROOT/'scripts/iso_irrep_exact.py',ROOT/'scripts/generate_irrep_data.py']}
print(json.dumps(dict(scope='SG5 four native scalar beta source models; bounded Q(zeta24) source certificate; no production or child-endpoint execution',specs=SPECS,counts=dict(counts),named_character_period=2,momentum_period=2,parent_collision_step=1,c_lattice='2x,2y,z integers; x-y integer',reciprocal='h,k,l integers; h+k even',input_sha256=pins),indent=2))
