import functools
from http.client import HTTPConnection
from http.server import ThreadingHTTPServer
import json
from pathlib import Path
import tempfile
import threading
import unittest

from serve_native_workbench import NativeWorkbenchHandler


class GatewayTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        (self.root / 'worker.json').write_text(json.dumps({'status': 'ready', 'limits': {'max_steps_per_command': 16}}))
        self.server = ThreadingHTTPServer(('127.0.0.1', 0), functools.partial(NativeWorkbenchHandler, directory=self.root))
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()

    def tearDown(self):
        self.server.shutdown()
        self.server.server_close()
        self.thread.join()
        self.temp.cleanup()

    def request(self, method='GET', path='/api/status', body=None, headers=None):
        connection = HTTPConnection(*self.server.server_address)
        payload = json.dumps(body) if body is not None else None
        supplied = {'Content-Type': 'application/json'} if payload is not None else {}
        supplied.update(headers or {})
        connection.request(method, path, payload, supplied)
        response = connection.getresponse()
        status, data = response.status, response.read()
        connection.close()
        return status, data

    def test_completion_ack_required_before_another_command(self):
        code, raw = self.request('POST', '/api/control', {'action': 'step', 'branch': 'both', 'steps': 4})
        self.assertEqual(code, 202)
        command = json.loads(raw)
        self.assertEqual(json.loads((self.root / 'command.json').read_text()), command)
        self.assertEqual(self.request('POST', '/api/control', {'action': 'observe'})[0], 409)
        (self.root / 'control.json').write_text(json.dumps({'id': command['id'], 'ok': True, 'status': 'completed'}))
        self.assertEqual(self.request('POST', '/api/control', {'action': 'observe'})[0], 202)
        self.assertEqual(len((self.root / 'gateway-requests.jsonl').read_text().splitlines()), 2)

    def test_bounds_and_closed_vocabulary(self):
        invalid = [{'action': 'eval', 'code': '1+1'}, {'action': 'step', 'steps': True},
                   {'action': 'step', 'steps': 17}, {'action': 'step', 'steps': 1.5},
                   {'action': 'observe', 'steps': 1}, {'action': 'step', 'steps': 1, 'branch': '../old'},
                   {'action': 'step', 'steps': 1, 'wall_seconds': 1000000}]
        for body in invalid:
            with self.subTest(body=body):
                self.assertEqual(self.request('POST', '/api/control', body)[0], 400)
        self.assertFalse((self.root / 'command.json').exists())

    def test_receptor_scale_contract_and_event_paths(self):
        for scales in ([0, 1, 1], [True, 1, 1], [1, 2], [float("inf"), 1, 1], [1e100, 1, 1], [1e-100, 1, 1]):
            self.assertEqual(self.request('POST', '/api/control', {'action': 'receptor_scales', 'branch': 'sensing', 'scales': scales})[0], 400)
        self.assertEqual(self.request('POST', '/api/control', {'action': 'receptor_scales', 'branch': 'off', 'scales': [1, 2, 3]})[0], 400)
        self.assertEqual(self.request('POST', '/api/control', {'action': 'observe', 'scales': [1, 2, 3]})[0], 400)
        self.assertEqual(self.request('POST', '/api/control', {'action': 'receptor_scales', 'branch': 'sensing', 'scales': [.5, .03, .006]})[0], 202)
        for phase in ('before', 'after'):
            frame = self.root / f'branches/sensing/frames/00001-{phase}'
            frame.mkdir(parents=True)
            (frame / 'frame.json').write_text('{}')
            self.assertEqual(self.request(path=f'/branches/sensing/frames/00001-{phase}/frame.json')[0], 200)
        self.assertEqual(self.request(path='/branches/sensing/frames/00001-untrusted/frame.json')[0], 404)

    def test_origin_and_host_checks(self):
        body = {'action': 'step', 'steps': 1}
        self.assertEqual(self.request('POST', '/api/control', body, {'Origin': 'https://example.invalid'})[0], 403)
        self.assertEqual(self.request(headers={'Host': 'example.invalid'})[0], 403)
        self.assertFalse((self.root / 'command.json').exists())

    def test_only_public_observations_are_served(self):
        frame = self.root / 'branches/off/frames/00001'
        frame.mkdir(parents=True)
        (frame / 'state.bin').write_bytes(b'actual-array')
        self.assertEqual(self.request(path='/branches/off/frames/00001/state.bin'), (200, b'actual-array'))
        for path in ['/worker.json', '/checkpoints/state.jls', '/branches/off/../../worker.json', '/branches/off/frames/00001/secret.txt']:
            self.assertEqual(self.request(path=path)[0], 404)
        (frame / 'frame.json').symlink_to(Path(__file__).resolve())
        self.assertEqual(self.request(path='/branches/off/frames/00001/frame.json')[0], 403)

    def test_starting_state_and_compact_embedded_status(self):
        (self.root / "worker.json").write_text(json.dumps({"status": "starting"}))
        self.assertEqual(json.loads(self.request()[1])["status"], "starting")
        self.assertEqual(self.request("POST", "/api/control", {"action": "describe"})[0], 409)
        status, html = self.request(path="/branches/off/index.html")
        self.assertEqual(status, 200)
        self.assertIn(b'id="workbench-embed"', html)
        self.assertIn(b'id="status"', html)
        self.assertIn(b'id="error"', html)

    def test_stopped_worker_and_journal_budget(self):
        (self.root / 'worker.json').write_text(json.dumps({'status': 'stopped'}))
        self.assertEqual(self.request('POST', '/api/control', {'action': 'describe'})[0], 409)
        (self.root / 'worker.json').write_text(json.dumps({'status': 'ready', 'limits': {'max_steps_per_command': 16}}))
        (self.root / 'gateway-requests.jsonl').write_text('x' * 2_000_000)
        self.assertEqual(self.request('POST', '/api/control', {'action': 'observe'})[0], 429)
        self.assertFalse((self.root / 'command.json').exists())


if __name__ == '__main__':
    unittest.main()
