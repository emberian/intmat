#!/usr/bin/env python3
"""Small HTTP integration checks for the external recorder."""
import gzip
import http.server
import json
from pathlib import Path
import subprocess
import tempfile
import threading
import unittest


class RecorderTest(unittest.TestCase):
    def test_http_observations(self):
        state = {'step': 1, 'time': .1, 'metrics': {'population': 4},
                 'agents': [{'genome_signature': 'a', 'energy': 2}, {'genome_signature': 'b', 'energy': 4}]}
        class Handler(http.server.BaseHTTPRequestHandler):
            def do_GET(self):
                if self.path == '/missing':
                    self.send_error(404)
                    return
                self.send_response(200)
                self.end_headers()
                self.wfile.write(json.dumps(state).encode())
            def log_message(self, *args):
                pass
        server = http.server.ThreadingHTTPServer(('127.0.0.1', 0), Handler)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            with tempfile.TemporaryDirectory() as directory:
                output = Path(directory) / 'run'
                base = f'http://127.0.0.1:{server.server_port}'
                cmd = ['python3', str(Path(__file__).with_name('record_arena_observations.py')),
                       '--source', f'ok={base}/latest.json', '--source', f'missing={base}/missing',
                       '--output', str(output), '--polls', '2', '--interval', '.01']
                subprocess.run(cmd, check=True)
                rows = [json.loads(row) for row in (output/'ok.jsonl').read_text().splitlines()]
                self.assertEqual(len(rows), 1)
                self.assertFalse(rows[0]['agents_complete'])
                self.assertEqual(rows[0]['captured_traits']['energy']['mean'], 3)
                self.assertEqual(len((output/'errors.jsonl').read_text().splitlines()), 2)
                archives = list(output.glob('*.gz'))
                self.assertEqual(len(archives), 1)
                self.assertEqual(json.loads(gzip.decompress(archives[0].read_bytes())), state)
                self.assertEqual(json.loads((output/'status.json').read_text())['state'], 'completed')
                refused = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
                self.assertNotEqual(refused.returncode, 0)
                self.assertEqual(len((output/'ok.jsonl').read_text().splitlines()), 1)
        finally:
            server.shutdown()
            server.server_close()


if __name__ == '__main__':
    unittest.main()
