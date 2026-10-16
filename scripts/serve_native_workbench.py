#!/usr/bin/env python3
"""Loopback gateway for an owned, bounded Julia native workbench session."""
import argparse
from functools import partial
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import math
from pathlib import Path
import re
import struct
import tempfile
import threading
import time
from urllib.parse import urlsplit, unquote
import uuid

REPO = Path(__file__).resolve().parent.parent
ACTIONS = {'describe', 'observe', 'step', 'pause', 'checkpoint', 'receptor_scales', 'material_rates', 'fork', 'stop'}
BRANCHES = {'off', 'sensing', 'fork1', 'fork2', 'both'}

def valid_scale(value):
    if type(value) not in (int, float):
        return False
    try:
        if not math.isfinite(value) or value <= 0:
            return False
        narrowed = struct.unpack('f', struct.pack('f', value))[0]
        return math.isfinite(narrowed) and narrowed > 0
    except (OverflowError, struct.error):
        return False

class NativeWorkbenchHandler(BaseHTTPRequestHandler):
    command_lock = threading.Lock()

    def __init__(self, *args, directory, **kwargs):
        self.root = Path(directory).resolve()
        super().__init__(*args, **kwargs)

    def log_request(self, code='-', size='-'):
        if self.command != 'GET':
            super().log_request(code, size)

    def local_request(self):
        host = self.headers.get('Host', '')
        return host in {f'127.0.0.1:{self.server.server_port}', f'localhost:{self.server.server_port}'}

    def send_bytes(self, status, data, mime='application/json'):
        self.send_response(status)
        self.send_header('Content-Type', mime)
        self.send_header('Content-Length', str(len(data)))
        self.send_header('Cache-Control', 'no-store')
        self.send_header('X-Content-Type-Options', 'nosniff')
        self.end_headers()
        self.wfile.write(data)

    def reply(self, status, value):
        self.send_bytes(status, json.dumps(value, allow_nan=False).encode())

    def read_json(self, name):
        return json.loads((self.root / name).read_text())

    def do_GET(self):
        if not self.local_request():
            return self.reply(403, {'error': 'Loopback host required'})
        path = unquote(urlsplit(self.path).path)
        try:
            if path == '/api/status':
                return self.reply(200, self.read_json('worker.json'))
            if path == '/api/ack':
                ack = self.root / 'control.json'
                return self.reply(200, self.read_json('control.json') if ack.exists() else {})
            if path in {'/', '/index.html'}:
                target = REPO / 'web/native_workbench.html'
            elif path == '/native_cell_compare.js':
                target = REPO / 'web/native_cell_compare.js'
            elif re.fullmatch(r'/branches/(off|sensing|fork1|fork2)/index.html', path):
                target = REPO / 'web/native_observatory.html'
            elif re.fullmatch(r'/branches/(off|sensing|fork1|fork2)/(observatory\.json|frames/[0-9]+(?:-(?:before|after))?/(frame\.json|state\.bin))', path):
                target = (self.root / path.lstrip('/')).resolve()
                if not target.is_relative_to(self.root):
                    return self.reply(403, {'error': 'Observation path leaves the session'})
            else:
                return self.reply(404, {'error': 'Unknown observation or page'})
            mime = {'.html': 'text/html; charset=utf-8', '.js': 'text/javascript; charset=utf-8', '.json': 'application/json', '.bin': 'application/octet-stream'}[target.suffix]
            data = target.read_bytes()
            if re.fullmatch(r'/branches/(off|sensing|fork1|fork2)/index.html', path):
                # Keep the native status/error nodes; only compress embedded chrome.
                style = b'<style id="workbench-embed">header{padding:8px 14px;justify-content:flex-end;min-height:38px}header>div:first-child{display:none}header .badge{padding:4px 9px;margin:0}main{padding:14px;gap:18px}</style>'
                data = data.replace(b'</style>', b'</style>' + style, 1)
            self.send_bytes(200, data, mime)
        except FileNotFoundError:
            self.reply(404, {'error': 'Not published yet'})
        except (OSError, ValueError):
            self.reply(503, {'error': 'Publication temporarily unavailable'})

    def do_POST(self):
        host = self.headers.get('Host', '')
        if not self.local_request() or self.headers.get('Origin', f'http://{host}') != f'http://{host}' or self.headers.get('Content-Type', '').split(';')[0] != 'application/json':
            return self.reply(403, {'error': 'Local same-origin JSON required'})
        if self.path != '/api/control':
            return self.reply(404, {'error': 'Unknown endpoint'})
        try:
            length = int(self.headers.get('Content-Length', '0'))
            if not 0 < length <= 2048:
                raise ValueError('Request size must be 1–2048 bytes')
            body = json.loads(self.rfile.read(length))
            if not isinstance(body, dict) or set(body) - {'action', 'branch', 'steps', 'scales', 'rates'}:
                raise ValueError('Only action, branch, steps, scales and rates are accepted')
            action, branch = body.get('action'), body.get('branch', 'both')
            if action not in ACTIONS or branch not in BRANCHES:
                raise ValueError('Unknown action or branch')
            worker = self.read_json('worker.json')
            if worker.get('status') != 'ready':
                return self.reply(409, {'error': 'Julia worker is not ready', 'status': worker.get('status')})
            if action not in worker.get('actions', ACTIONS) or (branch != 'both' and 'branches' in worker and branch not in worker['branches']):
                raise ValueError('Worker does not expose that action or branch')
            if action == 'fork' and branch == 'both':
                raise ValueError('Fork requires one existing branch')
            if action == 'step':
                limit = worker['limits']['max_steps_per_command']
                steps = body.get('steps')
                if type(steps) is not int or not 1 <= steps <= limit:
                    raise ValueError(f'Steps must be an integer in 1–{limit}')
            elif 'steps' in body:
                raise ValueError('Steps apply only to step')
            if action == 'receptor_scales':
                scales = body.get('scales')
                if branch == 'both' or (branch != 'sensing' and 'fork' not in worker.get('actions', [])) or not isinstance(scales, list) or len(scales) != 3 or not all(valid_scale(x) for x in scales):
                    raise ValueError('Three positive finite scales for sensing are required')
            elif 'scales' in body:
                raise ValueError('Scales apply only to receptor_scales')
            if action == 'material_rates':
                selected = worker.get('branches', {}).get(branch, {})
                controls = selected.get('material_controls', {})
                rates = body.get('rates')
                if branch == 'both' or not selected.get('parent') or not controls.get('enabled') or not isinstance(rates, dict) or set(rates) != {'decay_rate', 'desorption_rate', 'raw_diffusion'}:
                    raise ValueError('Material rates require one explicit fork and three rates')
                caps = {'decay_rate': controls['max_decay_rate'], 'desorption_rate': controls['max_desorption_rate'], 'raw_diffusion': controls['max_raw_diffusion']}
                if not all(type(v) in (int, float) and math.isfinite(v) and 0 <= v <= caps[k] and (v == 0 or valid_scale(v)) for k, v in rates.items()):
                    raise ValueError('Material rates outside branch bounds')
            elif 'rates' in body:
                raise ValueError('Rates apply only to material_rates')
            command = dict(body, branch=branch, id=str(uuid.uuid4()), requested_at=time.time())
            with self.command_lock:
                previous = self.root / 'command.json'
                if previous.exists():
                    ack = self.read_json('control.json') if (self.root / 'control.json').exists() else {}
                    if self.read_json('command.json').get('id') != ack.get('id'):
                        return self.reply(409, {'error': 'A command awaits Julia acknowledgement'})
                journal = self.root / 'gateway-requests.jsonl'
                line = json.dumps(command, allow_nan=False) + '\n'
                if journal.exists() and journal.stat().st_size + len(line.encode()) > 2_000_000:
                    return self.reply(429, {'error': 'Request journal budget reached; preserve this session'})
                with journal.open('a') as stream:
                    stream.write(line)
                with tempfile.NamedTemporaryFile('w', dir=self.root, prefix='.command-', delete=False) as stream:
                    json.dump(command, stream, allow_nan=False)
                    temporary = Path(stream.name)
                temporary.replace(previous)
            self.reply(202, command)
        except (ValueError, TypeError, KeyError):
            self.reply(400, {'error': 'Invalid command or worker contract'})
        except OSError:
            self.reply(503, {'error': 'Session files unavailable'})

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_directory', type=Path)
    parser.add_argument('--port', type=int, default=8807)
    args = parser.parse_args()
    root = args.run_directory.resolve()
    if not (root / 'worker.json').is_file():
        parser.error('Start the bounded Julia worker first')
    server = ThreadingHTTPServer(('127.0.0.1', args.port), partial(NativeWorkbenchHandler, directory=root))
    print(f'Native workbench: http://127.0.0.1:{server.server_port}/', flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()

if __name__ == '__main__':
    main()
