import io

import pytest
import pandas as pd
import numpy as np

import simsi_transfer.transfer as transfer


def test_calculate_average_probabilities():
    test1 = ['S(0.326)S(0.185)S(0.489)PPPRK', 'SS(1)SPPPRK', 'S(0.01)S(0.19)S(0.8)PPPRK', 'S(0.4)S(0.6)SPPPRK']
    assert transfer.calculate_average_probabilities(test1) == 'S(0.184)S(0.494)S(0.322)PPPRK'

    test2 = [np.nan, np.nan, np.nan, 'S(0.4)S(0.6)SPPPRK']
    assert transfer.calculate_average_probabilities(test2) == 'S(0.4)S(0.6)SPPPRK'

    test3 = [np.nan]
    assert np.isnan(transfer.calculate_average_probabilities(test3))
