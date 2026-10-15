#!/usr/bin/env python3
"""Serve the recorded gallery, including byte-range requests for video seeking."""
import argparse
import functools
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import os
from pathlib import Path
import re


class GalleryHandler(SimpleHTTPRequestHandler):
    def send_head(self):
        self.byte_range = None
        requested = self.headers.get('Range', '')
        path = self.translate_path(self.path)
        # Unknown/multiple ranges and conditional ranges fall back to a full
        # response. The standard handler retains redirects and cache handling.
        match = re.fullmatch(r'bytes=(\d*)-(\d*)', requested)
        if not match or self.headers.get('If-Range') or not os.path.isfile(path):
            return super().send_head()
        try:
            source = open(path, 'rb')
        except OSError:
            self.send_error(404, 'File not found')
            return None
        size = os.fstat(source.fileno()).st_size
        first, last = match.groups()
        if first:
            start = int(first)
            end = min(int(last), size - 1) if last else size - 1
        else:
            count = int(last) if last else 0
            start, end = max(0, size - count), size - 1
        if not size or start > end or start >= size:
            source.close()
            self.send_response(416)
            self.send_header('Content-Range', f'bytes */{size}')
            self.send_header('Content-Length', '0')
            self.end_headers()
            return None
        self.byte_range = (start, end)
        self.send_response(206)
        self.send_header('Content-Type', self.guess_type(path))
        self.send_header('Content-Range', f'bytes {start}-{end}/{size}')
        self.send_header('Content-Length', str(end - start + 1))
        self.send_header('Last-Modified', self.date_time_string(os.fstat(source.fileno()).st_mtime))
        self.end_headers()
        return source

    def end_headers(self):
        self.send_header('Accept-Ranges', 'bytes')
        super().end_headers()

    def copyfile(self, source, output):
        if self.byte_range is None:
            return super().copyfile(source, output)
        start, end = self.byte_range
        source.seek(start)
        remaining = end - start + 1
        while remaining:
            block = source.read(min(remaining, 64 * 1024))
            if not block:
                break
            output.write(block)
            remaining -= len(block)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, default=Path('_site'))
    parser.add_argument('--port', type=int, default=8770)
    args = parser.parse_args()
    handler = functools.partial(GalleryHandler, directory=str(args.directory.resolve()))
    ThreadingHTTPServer(('127.0.0.1', args.port), handler).serve_forever()
