import unittest
import numpy as np

from Transform_Data import Encoder


class MyTestCase(unittest.TestCase):
    def test_encode(self):
        x = ['a', 'a', 'a', 'b', 'b', 'c']
        y = [1, 1, 0, 1, 0, 0]
        encoder = Encoder()
        encoder.fit(x, y)

        self.assertEqual(2/3, encoder.encode('a'))
        self.assertEqual(1/2, encoder.encode('b'))
        self.assertEqual(0, encoder.encode('c'))
        self.assertEqual(0, encoder.encode('z'))

        encoded = encoder.transform(['a', 'b', 'a', 'b', 'c', 'z', np.nan])
        self.assertEqual([2/3, 1/2, 2/3, 1/2, 0, 0, 0], encoded)

        self.assertEqual(3, encoder.get_nr_classes())


if __name__ == '__main__':
    unittest.main()
