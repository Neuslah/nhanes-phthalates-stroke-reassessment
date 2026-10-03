from pathlib import Path
import csv,json,hashlib,math,sys,collections
sys.stdout.reconfigure(encoding='utf-8')
R=Path(r'E:\0317-NHANES_phthalates');C=R/'00_项目控制与AI工作流';EV=C/'02_任务执行证据'
W=EV/'20260930-PostAuthority-Rebuild-LineA-v0.2'/'02_TASK_B_TABLE1_DESCRIPTIVE_EVIDENCE'
FA=EV/'20260925-Statistical-Correction-PhaseA'/'resume_cholesterol_proxy';O=W/'03_OUTPUT'/'R1_FINAL';Q=W/'04_QC'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest().upper()
def rd(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def wc(p,r):
 assert not p.exists()
 with p.open('w',encoding='utf-8-sig',newline='') as f:
  c=csv.DictWriter(f,fieldnames=list(r[0]));c.writeheader();c.writerows(r)
before=rd(W/'00_INPUT_IDENTITY'/'R1_preservation_before.csv')
for r in before:assert sha(r['path'])==r['sha256'],r['path']
for filename,col in [('table1_R1_execution_QC.csv','pass'),('independent_Table1_QC.csv','pass'),('full_precision_serialization_QC.csv','exact_numeric_roundtrip')]:
 assert all(x[col]=='TRUE' for x in rd(Q/filename)),filename
flow=FA/'participant_flow_counts.csv';selection=FA/'included_excluded_selection.csv'
assert sha(flow)=='77B8196F2C8ED8E3D8468FCD5F56ECA41BFEF3C38480C722DC39C44DA46ED8D8'
assert sha(selection)=='18C3C596792A2B1BFAFE5D02A4C80DF4EBC0228E481ABA96581257F4C6F8B990'
f=[x for x in rd(flow) if x['group']=='OVERALL']
expected=[80312,44790,43849,43781,13634,13425,13329,10528,9366]
assert [int(x['n']) for x in f]==expected
flowrows=[]
for i,x in enumerate(f):
 n=int(x['n']);prev='' if i==0 else expected[i-1]
 flowrows.append({'stage_order':i+1,'saved_stage':x['stage'],'remaining_N_saved':n,'previous_remaining_N':prev,
  'adjacent_arithmetic_difference_N':'' if i==0 else prev-n,'derived_count_basis':'NOT_APPLICABLE' if i==0 else 'ARITHMETIC_ONLY',
  'branch_role':'COMMON_CYCLE_BRANCH_FROM_PRIMARY_CC' if i==8 else 'PRIMARY_FLOW',
  'reason_scope':'Stage-defined adjacent remaining-N subtraction; not mutually exclusive missing-reason accounting',
  'source_path':str(flow),'source_sha256':sha(flow)})
wc(O/'participant_flow_manuscript_source.csv',flowrows)
s=rd(selection);overall=[x for x in s if x['framework']=='MAIN_8_CYCLES' and x['group']=='OVERALL']
assert overall and all(int(x['included_total'])==10528 and int(x['excluded_total'])==2801 for x in overall)
assert 10528+2801==13329
variables=sorted({x['variable'] for x in overall})
assert set(variables)=={'stroke','RIDAGEYR','sex','race_eth_f','edu_f','marital_f','pir','bmi','smoking_f','hypertension_f','diabetes_f','hyperlipidemia_f'}
for x in overall:
 assert 0<=int(x['included_nonmissing'])<=10528 and 0<=int(x['excluded_nonmissing'])<=2801
 assert math.isfinite(float(x['absolute_SMD']))
for v in variables:assert {x['weighting'] for x in overall if x['variable']==v}=={'UNWEIGHTED','SURVEY_WEIGHTED'}
stroke=[x for x in overall if x['variable']=='stroke'];assert len(stroke)==2
availability={v:{'levels':sorted({x['level'] for x in overall if x['variable']==v}),
 'weighting_methods':sorted({x['weighting'] for x in overall if x['variable']==v}),
 'saved_rows':sum(x['variable']==v for x in overall),
 'included_observed_denominators':sorted({int(x['included_nonmissing']) for x in overall if x['variable']==v}),
 'excluded_observed_denominators':sorted({int(x['excluded_nonmissing']) for x in overall if x['variable']==v})} for v in variables}
selection_info={'source_path':str(selection),'bytes':selection.stat().st_size,'sha256':sha(selection),
 'authority_scope':'MAIN_8_CYCLES / OVERALL included-vs-excluded complete-case selection descriptive evidence only',
 'main_preCC_N':13329,'included_CC_N':10528,'excluded_N':2801,'saved_primary_overall_rows':len(overall),
 'weighted_unweighted_SMD_available':True,'stroke_prevalence_available':True,'variable_specific_observed_denominators_available':True,
 'availability':availability,'saved_stroke_records':stroke,
 'SMD_recomputed':False,'final_S3_rows_or_layout_selected':False,
 'limitations':'No proof of MAR/MNAR or absence of selection bias. Twelve saved variable classes only; no alcohol/exposure concentration selection comparison added. Other source domains/groups retain original provenance and are outside this primary-overall scoped registration.'}
sp=O/'selection_audit_evidence_closure.json';assert not sp.exists();sp.write_text(json.dumps(selection_info,ensure_ascii=False,indent=2),encoding='utf-8')
warnings=rd(O/'table1_warnings.csv');cnt=collections.Counter(x['message'] for x in warnings)
assert len(cnt)==38 and set(cnt.values())=={32}
wc(Q/'table1_warning_summary.csv',[{'message':k,'count':v,'classification':'LOCKED_DOMAIN_LONELY_PSU_NOTIFICATION','option':'survey.lonely.psu=adjust;survey.adjust.domain.lonely=TRUE'} for k,v in sorted(cnt.items())])
checks=[{'check':'remaining_N_chain_exact','pass':True,'detail':'80312→44790→43849→43781→13634→13425→13329→10528;9366 common branch'},
 {'check':'adjacent_differences_arithmetic_only','pass':True,'detail':'No participant-ID reason analysis or selection rerun'},
 {'check':'selection_primary_scope','pass':True,'detail':f'{len(overall)} saved rows;12 variables;both weighting methods'},
 {'check':'selection_denominators_finite_saved_SMD','pass':True,'detail':'Only saved fields inspected; no recomputation'},
 {'check':'protected_before_registration','pass':True,'detail':f'{len(before)} unique files unchanged'}]
wc(Q/'flow_selection_closure_QC.csv',checks)
files=[p for p in O.iterdir() if p.is_file()]+[W/'02_EXECUTION'/'05_corrected_table1_descriptive.R',W/'02_EXECUTION'/'06_independent_table1_QC.R',W/'02_EXECUTION'/'07_full_precision_export.R',Path(__file__),W/'AUTHOR_ADJUDICATION_TABLE1_CATEGORICAL_TEST_20260930.md']
wc(Q/'R1_output_identity.csv',[{'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)} for p in files])
report=f'''# Task B-R1 Table 1 / flow / selection evidence QC

TABLE1_DESCRIPTIVE_REBUILD = PASS
FLOW_EVIDENCE = CLOSED_FOR_MANUSCRIPT_REBUILD
S3_SELECTION_EVIDENCE = CLOSED_FOR_MANUSCRIPT_REBUILD
AUTHOR_ADJUDICATION_APPLIED = YES
CATEGORICAL_TEST_METHOD_UNIFORM = YES
HISTORICAL_P_REPRODUCTION_REQUIRED = NO

## Table 1 identity and execution

Frozen input SHA25643826EA86EAA72F598D35BE1B36FC3BEA282B1BB0C255A54E8A0F19025C3FE16; primary CC10528/stroke385/no stroke10143; Female5222/198 and Male5306/187. Exact Female/Male primary SEQN-set identity to Step0 PASS; no participants added/deleted. Complete-case indicator identical to saved precc AND complete.cases of locked M3 covariates. Factor/reference levels use the frozen object's actual metadata. Alcohol uses alcohol_harmonized_f; diabetes Yes/No labels only; structural proxy values unchanged.

Positive phthalate-weight base21970 is constructed BEFORE CC subsetting, ph_weight/8, SDMVPSU/SDMVSTRA,nest=TRUE; exact base membership/order/weights/strata/PSU identity to saved Phase-A base8 PASS. Primary represented244 PSU/120 strata,df124,zero lonely strata. The clinical flow count13634 is a later clinical domain and is not the full survey base.

Three survey::svyttest continuous tests and nine survey::svychisq(statistic="F") categorical tests only. Actual test statistic, numerator/denominator df, full-precision P and formatted P saved. Descriptive t df=123; adjusted F df retained exactly from installed survey4.5 (including noninteger multi-category df). No regression DOMAIN_DF override. No historical P comparison as QC criterion.

Continuous mean=svymean, SD=sqrt(svyvar). Independent direct formula uses normalized ph_weight/8 and n/(n-1) correction matching installed svyvar. Counts and weighted proportions independently verified. All raw/canonical/display source links, historical variable/category order and bold reference labels PASS. 461 independent QC checks, plus execution checks and exact numeric CSV roundtrip PASS. Full-precision CSVs use17 significant digits and round-trip identical to binary RDS; display means/SD/proportions1 decimal and P3 decimals/<0.001. No P reverse-engineering.

## Warning review and implementation detail

{len(warnings)} saved notifications are exactly38 stroke-domain lonely strata ×32 summary-function evaluations; no different warning. Primary/base designs have zero lonely strata. Locked survey.lonely.psu=adjust and survey.adjust.domain.lonely=TRUE are unchanged. Notifications preserved, not suppressed; no warning-count cutoff introduced. Finite means/SD/proportions/test statistics/P and independent values PASS.

Initial 03 guard assumed numeric diabetes codes; frozen field is a labelled factor. It stopped before tests. Actual-type-equivalent guard repair preserved original script/log/QC. Initial 04 warning wrapper then stopped on a locked-rule domain-lonely notification; 05 records those notifications with unchanged package behavior and fails on different warnings. Partial metadata/logs retained; no completed inferential model or Table 1 object was overwritten.

No direct glm/svyglm/lm/mice/BH/exposure-model call in the execution script; static call inventory saved. Installed survey::svyttest internally implements the authorized two-sample mean comparison through a Gaussian survey fit. This package-internal computation is confined to the three author-authorized descriptive tests; no analyst-defined association/regression model was fitted. This detail is disclosed instead of claiming the package executes no fitting internally.

## Flow / selection

Flow source3634 bytes/{sha(flow)}; nine verified remaining-N stages retained. New adjacent differences are arithmetic only, with the last1162 explicitly a common-cycle restriction branch. The30147 and2801 differences are not reinterpreted as mutually exclusive missing causes. Borderline301 is not split into primary flow.

Selection source205343 bytes/{sha(selection)}; primary preCC13329,CC10528,excluded2801; {len(overall)} original MAIN_8_CYCLES/OVERALL records cover12 variables and both weighting methods, with finite saved SMDs, saved stroke prevalence and variable-specific observed denominators. No SMD/statistical recomputation, participant selection, final S3 row choice or layout decision. No MAR/MNAR or no-selection-bias inference.

## Scope / hashes / readiness

Full execution source, R sessionInfo, package versions, installed method source, all outputs/bytes/SHA listed in R1_output_identity.csv. {len(before)} pre-existing protected files unchanged before authority write, including31 public-repository paths, all submitted artifacts, original Task B STOP evidence and controls. Old STOP/spec CONFLICT files remain historical provenance; author adjudication closes their boundary prospectively.

READY_FOR_MINIMAL_SCOPED_PA_REGISTRATION = YES
NEW_INFERENTIAL_ANALYSIS = NO / NO_NEW_ASSOCIATION_HYPOTHESIS;12_AUTHORIZED_TABLE1_DESCRIPTIVE_TESTS_ONLY
REGRESSION_MODEL_FIT = NO / NO_SEPARATE_ASSOCIATION_OR_REGRESSION_ANALYSIS
MI_RERUN = NO
FLOW_RECOMPUTED = NO
S3_SMD_RECOMPUTED = NO
MANUSCRIPT_MODIFIED = NO
SUPPLEMENT_MODIFIED = NO
FIGURES_MODIFIED = NO
STROBE_MODIFIED = NO
TASK_C_STARTED = NO
'''
rp=Q/'Task-B-R1-Table1-Flow-Selection-QC.md';assert not rp.exists();rp.write_text(report,encoding='utf-8')
print(json.dumps({'closure_QC':'PASS','selection_primary_overall_rows':len(overall),'warnings':len(warnings),'warning_strata':len(cnt),'ready_for_scoped_registration':True}))
