import unittest
import numpy as np

from CatEncoder import CatEncoder


class MyTestCase(unittest.TestCase):
    def test_encode(self):
        x = ['a', 'a', 'a', 'b', 'b', 'c']
        y = [1, 1, 0, 1, 0, 0]
        encoder = CatEncoder(feature="dummy")
        encoder.fit(x, y)

        self.assertEqual((2/3)/(2/3+1/2), encoder.encode_float('a'))
        self.assertEqual((1/2)/(2/3+1/2), encoder.encode_float('b'))
        self.assertEqual(0, encoder.encode_float('c'))
        self.assertEqual(0, encoder.encode_float('z'))

        encoded = encoder.transform(['a', 'b', 'a', 'b', 'c', 'z', np.nan])
        self.assertTrue(np.allclose(np.divide([2/3, 1/2, 2/3, 1/2, 0, 0, 0], (2/3+1/2)), encoded))

        self.assertEqual(3, encoder.get_nr_classes())


if __name__ == '__main__':
    unittest.main()
