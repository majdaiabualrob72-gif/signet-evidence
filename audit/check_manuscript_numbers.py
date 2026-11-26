"""GENESIS VI: compare Table 3 of the final manuscript with tables regenerated from archived outputs (make_tables.py)."""
import re, sys
md=open('experiments/simulation/tables_simulation.md').read()
tab={}
for s in re.split(r'^## ',md,flags=re.M)[1:]:
    name=s.split('\n')[0].strip(); rows={}
    for l in s.split('\n'):
        if l.startswith('| ') and not l.startswith('| Method') and not l.startswith('|---'):
            c=[x.strip() for x in l.strip().strip('|').split('|')]; rows[c[0]]=c
    tab[name]=rows
ms=open(sys.argv[1]).read()
t3=ms[ms.index('**Table 3.'):ms.index('**Table 4.')]
smap={'Reference':'reference','Null effect':'null effect','Biases doubled':'magnitudes doubled','Dense small biases':'dense small biases',
 '15% of signs reversed':'15% of signs wrong','Omitted bias':'omitted bias','No design opposing benefit':'no design opposing benefit','No bias':'no bias',
 'One-signed, null effect':'one-signed literature, null effect','One-signed, null effect, biases doubled':'one-signed literature, null effect, magnitudes doubled',
 'Stress: extreme outlier':'extreme outlier','Stress: publication bias':'publication bias','Stress: heterogeneous populations':'heterogeneous populations',
 'Stress: missing controls and variances':'missing controls and variances','Stress: contradictory evidence':'contradictory evidence'}
cols=['RE','BA-U','BA-S','NCO-cal','SIGNET protocol','SIGNET R3 k=1']
ok=bad=0
for l in t3.split('\n'):
    if not l.startswith('| ') or l.startswith('| Scenario') or l.startswith('| ---'): continue
    c=[x.strip() for x in l.strip().strip('|').split('|')]
    sc=smap[c[0]]
    for j,meth in enumerate(cols):
        cov,ben=re.match(r'([\d.]+) \(([\d.]+)\)',c[j+1]).groups()
        r=tab[sc][meth]
        if r[1]==cov and r[5]==ben: ok+=1
        else: bad+=1; print('MISMATCH',c[0],meth,'ms',cov,ben,'regen',r[1],r[5])
print(f'Table 3 cells matching regenerated outputs: {ok}/{ok+bad}')
