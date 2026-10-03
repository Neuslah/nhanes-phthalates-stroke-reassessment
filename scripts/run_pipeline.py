"""Reproduce current registered analyses in an empty external workspace.

Historical data/model/partial-output folders are never used as runtime inputs.
"""
from pathlib import Path
import argparse,csv,json,hashlib,subprocess,os,sys,shutil
sys.dont_write_bytecode=True
from compare_outputs import compare
from project_participant_flow import project
from membership_fingerprints import fingerprint

def rows(p):return list(csv.DictReader(Path(p).open(encoding='utf-8-sig')))
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest().upper()
def write_csv(p,z):
 with Path(p).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(z[0]));w.writeheader();w.writerows(z)
def write_json(p,z):Path(p).write_text(json.dumps(z,indent=2)+'\n',encoding='utf-8')

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--package-root',type=Path,default=Path(__file__).resolve().parents[1])
 parser.add_argument('--input-root',type=Path,required=True)
 parser.add_argument('--work-root',type=Path,required=True)
 parser.add_argument('--rscript',default='Rscript')
 a=parser.parse_args();package=a.package_root.resolve();inputs=a.input_root.resolve();work=a.work_root.resolve()
 if work==package or package in work.parents:raise RuntimeError('Participant-level workspace must be outside the public package')
 if work.exists() and any(work.iterdir()):raise RuntimeError('Refuse a nonempty workspace; this entrypoint is a fresh run')
 work.mkdir(parents=True,exist_ok=True)
 for n in ['resume_released_values','resume_cholesterol_proxy','analysis','identity','mi','logs','comparisons','tmp']:(work/n).mkdir()
 env=os.environ.copy()
 if os.name=='nt':
  for k in list(env):
   if k.startswith('LC_') or k in ['LANG','LANGUAGE']:env.pop(k)
  env['LANG']='English_United States.utf8'
 for k in ['TEMP','TMP','TMPDIR']:env[k]=str(work/'tmp')
 commands=[]
 def R(name,*args):
  cmd=[a.rscript,'--vanilla',str(package/'R'/name),*[str(p) for p in args]]
  log=work/'logs'/(name+'.log')
  commands.append(dict(step=name,command=cmd,cwd=str(work),log=str(log)))
  write_json(work/'commands.json',commands)
  print('RUN',name,flush=True)
  with log.open('wb') as f:result=subprocess.run(cmd,cwd=work,env=env,stdout=f,stderr=subprocess.STDOUT)
  if result.returncode:raise RuntimeError('R execution stopped: '+name+'; inspect '+str(log))
 def check(pa,name,actual,keys):
  compare(package/'expected'/pa/name,actual,work/'comparisons'/(pa+'_'+name),keys.split(','),'Registered '+pa+'; config/expected-authority-map.csv')
 try:
  # The registered-cache route has exact identities for every source file.
  alternatives={z['relative_path']:z['validated_fresh_acquisition_sha256'] for z in rows(package/'config/validated-acquisition-sha256.csv')}
  for z in rows(package/'config/source-inputs.csv'):
   p=inputs/z['relative_path']
   allowed={z['sha256'].upper()}
   if z['kind']=='translated_component':allowed.add(alternatives[z['relative_path']])
   if not p.is_file() or digest(p) not in allowed:raise RuntimeError('Source identity mismatch: '+z['relative_path'])
  protected=[dict(path=str(p),bytes=p.stat().st_size,sha256=digest(p)) for root in [inputs,package/'R'] for p in root.rglob('*') if p.is_file()]
  write_csv(work/'current_input_manifest.csv',protected)
  R('00_environment.R',package,work)
  R('01_source_and_alcohol.R',work/'resume_released_values',inputs)
  R('02_reconstruct_data.R',work/'resume_released_values',inputs,package)
  R('03_cholesterol_and_selection.R',work/'resume_cholesterol_proxy')
  R('03_export_core_identity.R',work/'resume_cholesterol_proxy',work/'identity')
  compare(package/'expected/core_counts.csv',work/'identity/core_counts.csv',work/'comparisons/core_counts.csv',['framework','sex'],'Registered core sample counts')
  fingerprint(work/'identity',work/'identity/domain_membership_sha256.csv')
  compare(package/'expected/domain_membership_sha256.csv',work/'identity/domain_membership_sha256.csv',work/'comparisons/domain_membership_sha256.csv',['domain'],'Registered Step0 domain membership fingerprints; no participant identifiers in package')
  project(work/'resume_cholesterol_proxy/participant_flow_counts.csv',work/'identity/participant_flow_manuscript_source.csv')
  check('PA-040','participant_flow_manuscript_source.csv',work/'identity/participant_flow_manuscript_source.csv','stage_order')
  check('PA-041','included_excluded_selection.csv',work/'resume_cholesterol_proxy/included_excluded_selection.csv','framework,group,variable,level,weighting')
  shutil.copyfile(package/'R/01_corrected_helpers.R',work/'analysis/01_corrected_helpers.R')
  R('04_primary.R',work/'analysis')
  check('PA-026','primary_all_models.csv',work/'analysis/results/primary_all_models.csv','analysis,sex,exposure,model')
  R('05_secondary_survey.R',work/'analysis')
  for pa,n,k in [('PA-028','common_cycle.csv','analysis,sex,exposure,model'),('PA-029','cycle_adjusted.csv','analysis,sex,exposure,model'),('PA-030','sex_interaction.csv','exposure'),('PA-031','RCS.csv','sex,exposure'),('PA-032','grouped.csv','sex,exposure'),('PA-033','male_subgroup_global.csv','exposure,modifier'),('PA-034','era_stratified.csv','sex,exposure,era'),('PA-034','era_interaction.csv','sex,exposure')]:check(pa,n,work/'analysis/results'/n,k)
  (work/'table1/03_OUTPUT/R1_FINAL').mkdir(parents=True);(work/'table1/04_QC').mkdir()
  R('08a_Table1.R',work/'table1',package);R('08b_Table1_serialize.R',work/'table1')
  for n in ['table1_canonical_full_precision.csv','table1_display.csv']:check('PA-039',n,work/'table1/03_OUTPUT/R1_FINAL'/n,'row_order')
  for n,k in [('table1_raw_descriptives_full_precision.csv','variable,level,group'),('table1_test_details_full_precision.csv','variable')]:check('PA-039',n,work/'table1/03_OUTPUT/R1_FINAL'/n,k)
  R('06a_MI_frameworks.R',work/'mi');R('06b_MI_initialize.R',work/'mi')
  if (work/'mi/precheck_blocking_issues.csv').exists() or (work/'mi/initialization_status.txt').read_text().strip()!='INITIALIZATION_PASS':raise RuntimeError('MI framework or initialization boundary')
  (work/'mi/completed').mkdir()
  for key in ['Female_eight_2003_2018','Male_eight_2003_2018','Female_ten_2005_2018','Male_ten_2005_2018']:
   R('06c_MI_FCS.R',work/'mi',key)
   if (work/'mi/completed'/key/'guarded_run_status.txt').read_text().strip()!='ALL_FCS_COMPUTED_PENDING_DIAGNOSTIC_REVIEW':raise RuntimeError('MI FCS stopped: '+key)
  R('06d_MI_identity.R',work/'mi',work/'identity');R('06e_MI_models.R',work/'analysis')
  check('PA-035','MI_Model3.csv',work/'analysis/results/MI_Model3.csv','analysis,sex,exposure,model')
  R('07_mixtures_no_RH.R',work/'analysis')
  for pa,n,k in [('PA-036','WQS_effects.csv','key'),('PA-036','WQS_weights.csv','sex,direction,component'),('PA-037','qgcomp_effects.csv','key'),('PA-037','qgcomp_component_weights.csv','sex,component')]:check(pa,n,work/'analysis/results'/n,k)
  write_json(work/'run_status.json',dict(status='PASS',RH_WQS_executed=False,new_analysis=False))
 except Exception as e:
  write_json(work/'run_status.json',dict(status='STOPPED',reason=str(e),public_release_ready=False))
  raise
if __name__=='__main__':main()
