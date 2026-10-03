from pathlib import Path
import csv,json,math,sys

ABS_TOL=1e-10
REL_TOL=1e-8
EXACT_NUMERIC={'N','events','DOMAIN_DF','overall_df','nonlinear_df','numerator_df','n','stroke_N','no_stroke_N','stroke_count','no_stroke_count','included_total','excluded_total','included_nonmissing','excluded_nonmissing','m','df_complete','rank','row_order'}

def records(p):return list(csv.DictReader(Path(p).open(encoding='utf-8-sig')))
def missing(v):return v in ('','NA','NaN','nan','None')
def compare(expected,actual,output,keys,source):
    a,b=records(expected),records(actual)
    checks=[]
    def add(row,col,ev,av,diff,tol,ok):checks.append(dict(row_key=row,column=col,expected_value=ev,clean_run_value=av,absolute_or_exact_difference=diff,tolerance=tol,status='PASS' if ok else 'FAIL',authority_source=source))
    if not a or not b:
        add('STRUCTURE','rows',len(a),len(b),abs(len(a)-len(b)),0,False)
    elif list(a[0])!=list(b[0]):
        add('STRUCTURE','columns',str(list(a[0])),str(list(b[0])),'EXACT_MISMATCH',0,False)
    else:
        add('STRUCTURE','row_count',len(a),len(b),abs(len(a)-len(b)),0,len(a)==len(b))
        am={tuple(x[k] for k in keys):x for x in a};bm={tuple(x[k] for k in keys):x for x in b}
        add('STRUCTURE','unique_keys',len(a),len(am),'EXACT',0,len(a)==len(am) and len(b)==len(bm))
        add('STRUCTURE','key_set',str(len(am)),str(len(bm)),'EXACT' if set(am)==set(bm) else 'EXACT_MISMATCH',0,set(am)==set(bm))
        for key in sorted(set(am)&set(bm)):
            for col,ev in am[key].items():
                av=bm[key][col];label='|'.join(key)
                if missing(ev) or missing(av):add(label,col,ev,av,'BOTH_MISSING' if missing(ev) and missing(av) else 'MISSING_MISMATCH',0,missing(ev) and missing(av));continue
                try:en,an=float(ev),float(av)
                except ValueError:add(label,col,ev,av,'EXACT' if ev==av else 'EXACT_MISMATCH',0,ev==av);continue
                exact=col in EXACT_NUMERIC or col in keys
                tolerance=0.0 if exact else ABS_TOL+REL_TOL*abs(en)
                diff=abs(en-an)
                add(label,col,ev,av,diff,tolerance,math.isfinite(en) and math.isfinite(an) and diff<=tolerance)
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(checks[0]));w.writeheader();w.writerows(checks)
    fails=[x for x in checks if x['status']=='FAIL']
    print(json.dumps(dict(expected=str(expected),actual=str(actual),checks=len(checks),failures=len(fails),first_failure=fails[:1]),ensure_ascii=False),flush=True)
    if fails:raise RuntimeError('STOP_NUMERIC_OR_IDENTITY_MISMATCH: '+str(output))
    return len(checks)

if __name__=='__main__':compare(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4].split(','),sys.argv[5])
