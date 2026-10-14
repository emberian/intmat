#!/usr/bin/env python3
"""Serve one arena run locally with bounded, auditable runtime controls."""
import argparse
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
import tempfile
import threading
import time
import uuid

ALLOWED = {'dt': (1e-7, 10), 'source_rate': (0, 100), 'sink_rate': (0, 100),
           'snapshot_seconds': (.1, 3600), 'checkpoint_seconds': (1, 86400)}

class ArenaHandler(SimpleHTTPRequestHandler):
    command_lock = threading.Lock()

    def end_headers(self):
        self.send_header('Cache-Control', 'no-store')
        super().end_headers()

    def reply(self, status, body):
        data = json.dumps(body).encode()
        self.send_response(status)
        self.send_header('Content-Type', 'application/json')
        self.send_header('Content-Length', str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        if self.path.split('?')[0] == '/api/status':
            self.reply(200, {'live_server': True, 'controls': list(ALLOWED)})
        else:
            super().do_GET()

    def do_POST(self):
        if self.path != '/api/control':
            return self.reply(404, {'error': 'Unknown endpoint'})
        # Requests must originate from this server; JSON prevents ordinary HTML form POSTs.
        host = self.headers.get('Host', '')
        origin = self.headers.get('Origin')
        if host not in {f'127.0.0.1:{self.server.server_port}', f'localhost:{self.server.server_port}'} or (origin and origin not in {f'http://{host}'}) or self.headers.get('Content-Type', '').split(';')[0] != 'application/json':
            return self.reply(403, {'error': 'Local same-origin JSON required'})
        try:
            length = int(self.headers.get('Content-Length', '0'))
            if not 0 < length <= 4096:
                raise ValueError('Invalid request size')
            body = json.loads(self.rfile.read(length))
            if not isinstance(body, dict):
                raise ValueError('Request must be an object')
            action = body.get('action')
            if action not in {'pause', 'resume', 'stop', 'checkpoint', 'set'}:
                raise ValueError('Unknown action')
            values = body.get('values', {})
            if not isinstance(values, dict) or (action != 'set' and values):
                raise ValueError('Values only apply to set')
            for key, value in values.items():
                if key not in ALLOWED or isinstance(value, bool) or not isinstance(value, (int, float)) or not ALLOWED[key][0] <= value <= ALLOWED[key][1]:
                    raise ValueError(f'Invalid value for {key}')
            command = {'id': str(uuid.uuid4()), 'action': action, 'values': values, 'requested_at': time.time()}
            root = Path(self.directory)
            with self.command_lock:
                command_file, ack_file = root / 'command.json', root / 'control.json'
                if command_file.exists():
                    previous = json.loads(command_file.read_text())
                    ack = json.loads(ack_file.read_text()) if ack_file.exists() else {}
                    if previous.get('id') != ack.get('id'):
                        return self.reply(409, {'error': 'A previous command is awaiting runtime acknowledgement'})
                with (root / 'control-requests.jsonl').open('a') as journal:
                    journal.write(json.dumps(command) + '\n')
                with tempfile.NamedTemporaryFile('w', dir=root, prefix='.command-', delete=False) as handle:
                    json.dump(command, handle)
                    temporary = Path(handle.name)
                temporary.replace(command_file)
            self.reply(202, command)
        except (ValueError, TypeError, json.JSONDecodeError) as error:
            self.reply(400, {'error': str(error)})

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_directory', type=Path)
    parser.add_argument('--port', type=int, default=8791)
    args = parser.parse_args()
    root = args.run_directory.resolve()
    if not root.is_dir():
        parser.error('run directory does not exist')
    from build_arena_view import build
    build(root)
    server = ThreadingHTTPServer(('127.0.0.1', args.port), partial(ArenaHandler, directory=str(root)))
    print(f'Arena: http://127.0.0.1:{args.port}/arena.html', flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        server.server_close()

if __name__ == '__main__':
    main()
