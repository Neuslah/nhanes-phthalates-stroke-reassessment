from pathlib import Path
from decimal import Decimal
import csv,hashlib,sys
def fingerprint(root,dest):
 z=[]
 for p in sorted(Path(root).glob('*_SEQN.txt')):
  values=[Decimal(x) for x in p.read_text(encoding='utf-8-sig').split()]
  if any(x!=x.to_integral_value() for x in values):raise RuntimeError('Noninteger participant identifier')
  values=sorted(int(x) for x in values)
  if len(values)!=len(set(values)):raise RuntimeError('Duplicate participant identifier')
  payload=''.join(str(x)+'\n' for x in values).encode('utf-8')
  z.append(dict(domain=p.stem.removesuffix('_SEQN'),N=len(values),canonical_SEQN_sha256=hashlib.sha256(payload).hexdigest().upper()))
 with Path(dest).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(z[0]));w.writeheader();w.writerows(z)
if __name__=='__main__':fingerprint(sys.argv[1],sys.argv[2])
