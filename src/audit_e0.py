"""Run E0 and bind its acceptance to the tested numerical sources."""
import hashlib,json,sys,unittest
from pathlib import Path
from model import Parameters,simulate,weights
ROOT=Path(__file__).resolve().parents[1]

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    suite=unittest.defaultTestLoader.discover(str(ROOT/'tests'))
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    old=json.loads((ROOT/'qa/E0_report.json').read_text())
    benchmark=json.loads((ROOT/'qa/benchmark_szabo.json').read_text())
    benchmark_unchanged=old['source_sha256']['src\\benchmark.py']==digest(ROOT/'src/benchmark.py')
    if not benchmark_unchanged:raise RuntimeError('Rerun canonical benchmark: source changed')
    checks=[]
    for eps in [0.,.05]:
        for kernel in ['W','E','PP']:
            p=Parameters(kernel=kernel,L=50,b=1.05,M=65,a=.4,m=8,epsilon=eps)
            r=simulate(p,burn=0,measure=100,every=5,seed=58,early=False)['row']
            checks.append(dict(kernel=kernel,epsilon=eps,L=50,G=r['G'],G_from_Cy=r['G_from_Cy'],
                relative_error=r['identity_absolute_error']/abs(r['G']),tail_weight=r['kernel_tail_weight'],
                tolerance=r['identity_tolerance'],n_pairs=r['n_pairs'],maximum_lag=len(weights(p)[0])-1))
    passed=result.wasSuccessful() and benchmark['passed'] and benchmark_unchanged
    sources=[ROOT/'src/model.py',ROOT/'src/benchmark.py',ROOT/'tests/test_model.py',ROOT/'config/E0.json']
    report=dict(E0_accepted=passed,implementation_tests=dict(passed=result.testsRun,failures=len(result.failures),errors=len(result.errors)),
        identity_checks=checks,canonical_benchmark=benchmark,canonical_benchmark_reused='New simple-paper benchmark already completed; source unchanged',
        benchmark_csv_sha256=digest(ROOT/'qa/benchmark_szabo.csv'),
        source_sha256={str(p.relative_to(ROOT)):digest(p) for p in sources})
    (ROOT/'qa/E0_report.json').write_text(json.dumps(report,indent=2))
    (ROOT/'qa/E0_STATUS.txt').write_text('E0 PASSED; exact observed-record identity, pooled coexistence sample.\n' if passed else 'E0 FAILED\n')
    gate=ROOT/'config/approvals.json';state=json.loads(gate.read_text());state['E0_passed']=passed
    gate.write_text(json.dumps(state,indent=2))
    print(json.dumps(report,indent=2))
    if not passed:sys.exit(1)
if __name__=='__main__':main()
