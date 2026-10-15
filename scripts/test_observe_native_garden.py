#!/usr/bin/env python3
"""Observer admission/failure tests; byte fixtures never become displayed frames."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from observe_native_garden import Observer, completed_entries

class AdmissionTests(unittest.TestCase):
    def setUp(self):
        self.temporary=tempfile.TemporaryDirectory();self.root=Path(self.temporary.name)
        self.source=self.root/'run';(self.source/'frames').mkdir(parents=True)
        self.payload=b'complete fixture';(self.source/'frames/frame.jls').write_bytes(self.payload)
        self.entry=dict(file='frames/frame.jls',bytes=len(self.payload),sha256=hashlib.sha256(self.payload).hexdigest(),step=1,time=.1,side=8)
    def tearDown(self):self.temporary.cleanup()
    def test_admission_rejects_wrong_hash_size_and_path(self):
        self.assertEqual(len(completed_entries(self.source,{'frames':[self.entry]})),1)
        for change in [dict(sha256='0'*64),dict(bytes=0),dict(file='../frame.jls'),dict(file='/tmp/frame.jls')]:
            self.assertEqual(completed_entries(self.source,{'frames':[dict(self.entry,**change)]}),[])
    def test_partial_manifest_retains_last_good(self):
        observer=Observer(self.source,self.root/'view');observer.index['frames']=[{'source_sha256':'kept'}]
        (self.source/'run.json').write_text('{"frames":[')
        self.assertFalse(observer.scan());self.assertEqual(observer.index['frames'],[{'source_sha256':'kept'}])
        self.assertIn('retained',json.loads((self.root/'view/observatory.json').read_text())['observer']['status'])
    def test_budget_stops_before_export(self):
        (self.source/'run.json').write_text(json.dumps({'frames':[self.entry]}))
        observer=Observer(self.source,self.root/'view',budget=1,julia='/nonexistent')
        observer.scan();self.assertEqual(observer.index['frames'],[])
        self.assertEqual(observer.index['observer']['stop_reason'],'observer artifact budget')
    def test_same_basename_different_path_rejected(self):
        observer=Observer(self.source,self.root/'view')
        from observe_native_garden import atomic_json
        atomic_json(self.root/'view/observatory.json',observer.index)
        other=self.root/'other/run';other.mkdir(parents=True)
        with self.assertRaisesRegex(ValueError,'source path'):
            Observer(other,self.root/'view')
    def test_restart_archives_stop_and_rechecks_budget(self):
        (self.source/'run.json').write_text(json.dumps({'frames':[self.entry]}))
        observer=Observer(self.source,self.root/'view',budget=1);observer.scan()
        resumed=Observer(self.source,self.root/'view',budget=1000000,julia='/nonexistent')
        self.assertNotIn('stop_reason',resumed.index['observer'])
        self.assertEqual(resumed.index['observer']['previous_session_stop_reason'],'observer artifact budget')
        resumed.scan()
        self.assertIn('export failed',resumed.index['observer']['stop_reason'])
    def test_replaced_run_at_same_path_rejected(self):
        (self.source/'run.json').write_text(json.dumps({'frames':[self.entry]}))
        observer=Observer(self.source,self.root/'view',budget=1);observer.scan()
        (self.source/'run.json').write_text(json.dumps({'frames':[dict(self.entry,sha256='f'*64)]}))
        with self.assertRaisesRegex(ValueError,'first-frame'):
            Observer(self.source,self.root/'view').scan()
    def test_export_failure_retains_observations(self):
        (self.source/'run.json').write_text(json.dumps({'frames':[self.entry]}))
        observer=Observer(self.source,self.root/'view',julia='/nonexistent');observer.index['frames']=[{'source_sha256':'kept','time':0,'side':8}]
        observer.scan();self.assertEqual(len(observer.index['frames']),1)
        self.assertIn('export failed',observer.index['observer']['stop_reason'])
if __name__=='__main__':unittest.main()
