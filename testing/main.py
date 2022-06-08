import json
from rdkit import Chem

def rdkit_pattern_fingerprint_test(request):
    try:
        return_value = []
        request_json = request.get_json()
        calls = request_json['calls']
        
        for call in calls:     
            smiles = call[0]  
            try:
                fp = Chem.PatternFingerprint(Chem.MolFromSmiles(smiles), tautomerFingerprints=True )
                return_value.append(fp.ToBinary())
            except:
                return_value.append("")

        return_json = json.dumps( { "replies" :  return_value} ), 200
        return return_json
    except Exception:
        return json.dumps( { "errorMessage": 'something unexpected in input' } ), 400


