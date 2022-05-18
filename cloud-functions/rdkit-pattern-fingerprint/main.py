import json
import re
from rdkit import Chem

def rdkit_pattern_fingerprint(request):

    """Responds to any HTTP request.
    Args:
        request (flask.Request): HTTP request object.
    Returns:
        The response text or any set of values that can be turned into a
        Response object using
        `make_response <http://flask.pocoo.org/docs/1.0/api/#flask.Flask.make_response>`.
    """

    try:
        return_value = []
        request_json = request.get_json()
        calls = request_json['calls']
        
        for call in calls:
        
            smiles = call[0]
            
            try:
                fp = Chem.RDKFingerprint(Chem.MolFromSmiles(smiles))
                return_value.append(fp.ToBase64())
            except:
                return_value.append("")

        return_json = json.dumps( { "replies" :  return_value} ), 200
        return return_json
    except Exception:
        return json.dumps( { "errorMessage": 'something unexpected in input' } ), 400