from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdqueries
import json

def rdkit_pattern_fingerprint_test(request):
    try:
        return_value = []
        request_json = request.get_json()
        calls = request_json['calls']
        
        for call in calls:     
            
            smiles = call[0]  
            
            try:
                mol = Chem.MolFromSmiles(smiles)
                fp_pattern_short = Chem.PatternFingerprint(mol, 11, tautomerFingerprints=True)
                fp_pattern_short_as_binary = DataStructs.BitVectToBinaryText(fp_pattern_short)
                fp_pattern_short_as_binary_hex = fp_pattern_short_as_binary.hex()
                fp_pattern_short_as_int = int(fp_pattern_short.ToBitString(),2)

                fp_pattern_long = Chem.PatternFingerprint(mol, tautomerFingerprints=True)
                fp_pattern_long_as_binary = DataStructs.BitVectToBinaryText(fp_pattern_long)
                fp_pattern_long_as_binary_hex = fp_pattern_long_as_binary.hex()
                fp_pattern_long_as_int = int(fp_pattern_long.ToBitString(),2)

                fp_morgan = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius = 2, useChirality=True, nBits = 2048)
                fp_morgan_as_binary = DataStructs.BitVectToBinaryText(fp_morgan)
                fp_morgan_as_binary_hex = fp_morgan_as_binary.hex()
                fp_morgan_as_int = int(fp_morgan.ToBitString(),2)

                fingerprints = {}

                fingerprints["fp_pattern_short_as_binary_hex"] = fp_pattern_short_as_binary_hex
                fingerprints["fp_pattern_long_as_binary_hex"] = fp_pattern_long_as_binary_hex
                fingerprints["fp_morgan_as_binary_hex"] = fp_morgan_as_binary_hex

                fingerprints["fp_pattern_short_as_int"]= fp_pattern_short_as_int
                fingerprints["fp_pattern_long_as_int"] = fp_pattern_long_as_int
                fingerprints["fp_morgan_as_int"] = fp_morgan_as_int

                # count carbons
                q = rdqueries.AtomNumEqualsQueryAtom(6)
                num_carbon = len(mol.GetAtomsMatchingQuery(q))

                # count nitrogens
                q = rdqueries.AtomNumEqualsQueryAtom(7)
                num_nitrogen = len(mol.GetAtomsMatchingQuery(q))

                # count oxygens
                q = rdqueries.AtomNumEqualsQueryAtom(8)
                num_oxygen = len(mol.GetAtomsMatchingQuery(q))

                # count fluorines
                q = rdqueries.AtomNumEqualsQueryAtom(9)
                num_fluorine = len(mol.GetAtomsMatchingQuery(q))

                # count sulfurs
                q = rdqueries.AtomNumEqualsQueryAtom(16)
                num_sulfur = len(mol.GetAtomsMatchingQuery(q))

                fingerprints["num_carbon"] = num_carbon
                fingerprints["num_nitrogen"] = num_nitrogen
                fingerprints["num_oxygen"] = num_oxygen
                fingerprints["num_fluorine"] = num_fluorine
                fingerprints["num_sulfur"] = num_sulfur


                fingerprints_as_json_string = json.dumps(fingerprints)
                return_value.append(fingerprints_as_json_string)
            
            except:
                return_value.append("")

        return_json = json.dumps( { "replies" :  return_value} ), 200
        return return_json
    except Exception:
        return json.dumps( { "errorMessage": 'something unexpected in input' } ), 400


