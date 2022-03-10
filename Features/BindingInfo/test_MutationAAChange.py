from unittest import TestCase

from Features.BindingInfo.MutationAAChange import *


class TestMutationAAChange(TestCase):

    def test_constructor(self):
        mutationAAChange = MutationAAChange()

        self.assertTrue(mutationAAChange.get_aa_coeff('W', True) > mutationAAChange.get_aa_coeff('F', True) )
        self.assertTrue(mutationAAChange.get_aa_coeff('F', True) > mutationAAChange.get_aa_coeff('Y', True) )
        self.assertTrue(mutationAAChange.get_aa_coeff('D', True) > mutationAAChange.get_aa_coeff('H', True) )
        self.assertTrue(mutationAAChange.get_aa_coeff('Q', True) > mutationAAChange.get_aa_coeff('E', True) )

        self.assertEqual(mutationAAChange.get_aa_coeff('W', True),
                         mutationAAChange.binding_pos_fact*mutationAAChange.aa_prime_regr_coeff['W'])
