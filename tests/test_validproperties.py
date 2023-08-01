# -*- coding: utf-8 -*-

"""
Created on Tue 29 Jun 2021 19:00:00

@author: Mark.Illingworth
"""

import pytest

from Refactor.Properties import ValidProperty

valid_testlist = [
    (3.5, 1.0, 10.0, 3.5, True),  # conventional value in range
    (3.5, 10.0, 1.0, 3.5, True),  # limits swapped should still work
    ]
invalid_testlist = [
    (10.9, 1.0, 10.0, 10.0, False),  # outside upper limit
    (0.9, 10.0, 1.0, 1.0, False),      # outside upper limit
    ]


@pytest.mark.parametrize("value," +
                         "limit1," +
                         "limit2," +
                         "expected_value," +
                         "expected_valid", valid_testlist)
def test_validated_value_returns_expected(value, limit1, limit2,
                                          expected_value, expected_valid):
    vv = ValidProperty(value, limit1, limit2)
    assert expected_value == vv.value
    assert expected_valid == vv.InRange


@pytest.mark.parametrize("value," +
                         "limit1," +
                         "limit2," +
                         "expected_value," +
                         "expected_valid", invalid_testlist)
def test_invalid_value_raises_exception(value, limit1, limit2,
                                        expected_value, expected_valid):
    with pytest.raises(ValueError) as error:
        vv = ValidProperty(value, limit1, limit2)
    assert str(error.value) == "Value out of range"
