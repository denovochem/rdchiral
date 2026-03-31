import json
import os

import pytest

from rdchiral.template_extractor import extract_from_reaction


def _load_extraction_cases():
    with open(
        os.path.join(os.path.dirname(__file__), "test_extraction_cases.json"), "r"
    ) as fid:
        return json.load(fid)


_EXTRACTION_CASES = _load_extraction_cases()


@pytest.mark.slow
@pytest.mark.parametrize(
    "test_case",
    [
        pytest.param(
            case,
            marks=pytest.mark.xfail(
                reason="Known invalid extraction case",
                strict=False,
            ),
        )
        if str(case.get("valid", "True")).lower() != "true"
        else case
        for case in _EXTRACTION_CASES
    ],
    ids=[str(case.get("_id", f"case_{i}")) for i, case in enumerate(_EXTRACTION_CASES)],
)
def test_template_extraction_case(test_case):
    extracted = extract_from_reaction(dict(test_case))["reaction_smarts"]
    assert extracted == test_case["expected_template"]
