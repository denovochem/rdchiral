import json
import os
import sys
import time

sys.path = [os.path.dirname(os.path.dirname((__file__)))] + sys.path

from rdchiral.main import (
    rdchiralReactants,
    rdchiralReaction,
    rdchiralRun,
    rdchiralRunText,
)

with open(
    os.path.join(os.path.dirname(__file__), "test_rdchiral_cases.json"), "r"
) as fid:
    test_cases = json.load(fid)

start = time.time()
all_passed = True
for i, test_case in enumerate(test_cases):
    print("\n# Test {:2d}/{}".format(i + 1, len(test_cases)))

    # Directly use SMILES/SMARTS
    reaction_smarts = test_case["smarts"]
    reactant_smiles = test_case["smiles"]
    if rdchiralRunText(reaction_smarts, reactant_smiles) == test_case["expected"]:
        print("    from text: passed")
    else:
        print("    from text: failed")
        all_passed = False

    # Pre-initialize & repeat
    rxn = rdchiralReaction(reaction_smarts)
    reactants = rdchiralReactants(reactant_smiles)
    if all(rdchiralRun(rxn, reactants) == test_case["expected"] for j in range(3)):
        print("    from init: passed")
    else:
        print("    from init: failed")
        all_passed = False

all_passed = "All passed!" if all_passed else "Failed!"
print("\n# Final result: {}".format(all_passed))
print("Total time: {} s".format(time.time() - start))
