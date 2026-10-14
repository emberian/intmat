"""Checks the actual local command gateway, without running a simulation."""
import functools
import http.client
from http.server import ThreadingHTTPServer
import json
from pathlib import Path
import tempfile
import threading
import unittest
from serve_arena import ArenaHandler

class ControlsTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.server = ThreadingHTTPServer(('127.0.0.1', 0), functools.partial(ArenaHandler, directory=self.directory.name))
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()
    def tearDown(self):
        self.server.shutdown()
        self.server.server_close()
        self.directory.cleanup()
    def post(self, body, origin=None):
        connection = http.client.HTTPConnection('127.0.0.1', self.server.server_port)
        headers = {'Content-Type': 'application/json'}
        if origin:
            headers['Origin'] = origin
        connection.request('POST', '/api/control', json.dumps(body), headers)
        response = connection.getresponse()
        result = response.status, json.loads(response.read())
        connection.close()
        return result
    def test_atomic_command(self):
        status, result = self.post({'action': 'set', 'values': {'source_rate': 0.25}})
        self.assertEqual(status, 202)
        saved = json.loads((Path(self.directory.name) / 'command.json').read_text())
        self.assertEqual(saved, result)
        self.assertEqual(saved['values'], {'source_rate': 0.25})
    def test_journal_rotation(self):
        root = Path(self.directory.name)
        journal = root / 'control-requests.jsonl'
        journal.write_text('x' * (1024 * 1024))
        self.assertEqual(self.post({'action': 'pause'})[0], 202)
        self.assertEqual((root / 'control-requests.previous.jsonl').stat().st_size, 1024 * 1024)
        self.assertLess(journal.stat().st_size, 4096)

    def test_pending_command_is_not_overwritten(self):
        status, first = self.post({'action': 'pause'})
        self.assertEqual(status, 202)
        self.assertEqual(self.post({'action': 'resume'})[0], 409)
        self.assertEqual(json.loads((Path(self.directory.name) / 'command.json').read_text())['id'], first['id'])
        (Path(self.directory.name) / 'control.json').write_text(json.dumps({'id': first['id'], 'ok': True}))
        self.assertEqual(self.post({'action': 'resume'})[0], 202)

    def test_nonobject_rejected(self):
        self.assertEqual(self.post([])[0], 400)

    def test_cross_origin_rejected(self):
        self.assertEqual(self.post({'action': 'pause'}, 'https://outside.example')[0], 403)
        self.assertFalse((Path(self.directory.name) / 'command.json').exists())
    def test_invalid_parameters_rejected(self):
        for values in [{'dt': -1}, {'unknown': 2}, {'sink_rate': True}, {'dt': float('nan')}]:
            self.assertEqual(self.post({'action': 'set', 'values': values})[0], 400)
        self.assertFalse((Path(self.directory.name) / 'command.json').exists())

if __name__ == '__main__':
    unittest.main()
