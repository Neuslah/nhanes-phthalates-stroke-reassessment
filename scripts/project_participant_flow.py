"""Project registered participant-flow stages from newly reconstructed counts."""
from pathlib import Path
import csv,hashlib,sys
def project(source,dest):
 source=Path(source);z=list(csv.DictReader(source.open(encoding='utf-8-sig')))
 overall={x['stage']:int(x['n']) for x in z if x['group']=='OVERALL'}
 stages=['pooled','adult','nonpregnant','valid_stroke','positive_weight_design','creatinine_available','eight_analytes_available_preCC','covariates_complete','common_cycle_complete']
 result=[];previous=None
 for order,stage in enumerate(stages,1):
  n=overall[stage]
  result.append(dict(stage_order=order,saved_stage=stage,remaining_N_saved=n,previous_remaining_N='' if previous is None else previous,adjacent_arithmetic_difference_N='' if previous is None else previous-n,derived_count_basis='NOT_APPLICABLE' if previous is None else 'ARITHMETIC_ONLY',branch_role='COMMON_CYCLE_BRANCH_FROM_PRIMARY_CC' if order==9 else 'PRIMARY_FLOW',reason_scope='Stage-defined adjacent remaining-N subtraction; not mutually exclusive missing-reason accounting',source_path='resume_cholesterol_proxy/participant_flow_counts.csv',source_sha256=hashlib.sha256(source.read_bytes()).hexdigest().upper()))
  previous=n
 with Path(dest).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(result[0]));w.writeheader();w.writerows(result)
if __name__=='__main__':project(sys.argv[1],sys.argv[2])
