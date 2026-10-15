import functools
from http.client import HTTPConnection
from http.server import ThreadingHTTPServer
from pathlib import Path
import tempfile
import threading
import unittest

from serve_gallery import GalleryHandler


class QuietHandler(GalleryHandler):
    def log_message(self, *args):
        pass


class GalleryRanges(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        Path(self.temporary.name, 'film.mp4').write_bytes(bytes(range(256)))
        handler = functools.partial(QuietHandler, directory=self.temporary.name)
        self.server = ThreadingHTTPServer(('127.0.0.1', 0), handler)
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()

    def tearDown(self):
        self.server.shutdown()
        self.server.server_close()
        self.thread.join()
        self.temporary.cleanup()

    def request(self, range=None, method='GET', extra=None):
        connection = HTTPConnection(*self.server.server_address)
        headers = dict(extra or {})
        if range is not None:
            headers['Range'] = range
        connection.request(method, '/film.mp4', headers=headers)
        response = connection.getresponse()
        result = response.status, dict(response.getheaders()), response.read()
        connection.close()
        return result

    def test_closed_open_suffix_ranges(self):
        for requested, start, end in [('bytes=17-25', 17, 25), ('bytes=250-', 250, 255), ('bytes=-7', 249, 255), ('bytes=250-999', 250, 255)]:
            with self.subTest(requested=requested):
                status, headers, body = self.request(requested)
                self.assertEqual(status, 206)
                self.assertEqual(headers['Content-Range'], f'bytes {start}-{end}/256')
                self.assertEqual(body, bytes(range(start, end + 1)))

    def test_unsatisfiable_range(self):
        for requested in ('bytes=256-', 'bytes=90-20', 'bytes=-0', 'bytes=-'):
            status, headers, body = self.request(requested)
            self.assertEqual((status, headers['Content-Range'], body), (416, 'bytes */256', b''))

    def test_full_conditional_and_head(self):
        status, headers, body = self.request()
        self.assertEqual((status, body), (200, bytes(range(256))))
        self.assertEqual(headers['Accept-Ranges'], 'bytes')
        self.assertEqual(self.request('bytes=2-4', extra={'If-Range': 'unknown'})[0], 200)
        status, headers, body = self.request('bytes=2-4', method='HEAD')
        self.assertEqual((status, headers['Content-Length'], body), (206, '3', b''))


if __name__ == '__main__':
    unittest.main()
