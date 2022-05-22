import json
from rdkit import Chem

def rdkit_pattern_fingerprint(request):
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

def rdkit_draw_svg(request):
    from rdkit.Chem.Draw import rdMolDraw2D
    import re
    pattern = re.compile("<\?xml.*\?>")
    try:
        return_value = []
        request_json = request.get_json()
        calls = request_json['calls']
        
        for call in calls:       
            smiles = call[0]            
            try:
                mc = Chem.MolFromSmiles(smiles)
                drawer = rdMolDraw2D.MolDraw2DSVG(*(450,150))
                drawer.DrawMolecule(mc)
                drawer.FinishDrawing()
                svg = drawer.GetDrawingText().replace('svg:', '')
                svg = re.sub(pattern, '', svg)
                return_value.append(svg)
            except:
                return_value.append("")

        return_json = json.dumps( { "replies" :  return_value} ), 200
        return return_json
    except Exception:
        return json.dumps( { "errorMessage": 'something unexpected in input' } ), 400

def rdkit_generate_MACCS_keys(request):
    from rdkit.Chem import MACCSkeys
    try:
        return_value = []
        request_json = request.get_json()
        calls = request_json['calls']
        
        for call in calls:       
            smiles = call[0]            
            try:
                mol = Chem.MolFromSmiles(smiles)
                fp = MACCSkeys.GenMACCSKeys(mol)
                return_value.append(fp.ToBase64())
            except:
                return_value.append("")

        return_json = json.dumps( { "replies" :  return_value} ), 200
        return return_json
    except Exception:
        return json.dumps( { "errorMessage": 'something unexpected in input' } ), 400

def rdkit_smiles_to_inchi(request):
    try:
        return_value = []
        request_json = request.get_json()
        calls = request_json['calls']
        
        for call in calls:     
            smiles = call[0]  
            try:
                mol = Chem.MolFromSmiles(smiles)
                return_value.append(Chem.MolToInchi(mol))
            except:
                return_value.append("")

        return_json = json.dumps( { "replies" :  return_value} ), 200
        return return_json
    except Exception:
        return json.dumps( { "errorMessage": 'something unexpected in input' } ), 400

def rdkit_descriptor_generic(request):
    from rdkit.Chem import Descriptors
    try:
        return_value = []
        request_json = request.get_json()
        calls = request_json['calls']
        rdkit_func = request_json["userDefinedContext"]["rdkit-function"]
        #print(rdkit_func)
        for call in calls:     
            smiles = call[0]  
            try:
                mol = Chem.MolFromSmiles(smiles)
                #print(getattr(Descriptors, rdkit_func)(mol))
                return_value.append(str(getattr(Descriptors, rdkit_func)(mol)))
            except:
                return_value.append("")

        return_json = json.dumps( { "replies" :  return_value} ), 200
        return return_json
    except Exception:
        return json.dumps( { "errorMessage": 'something unexpected in input' } ), 400

def rdkit_qed(request):
    from qed import qed
    try:
        return_value = []
        request_json = request.get_json()
        calls = request_json['calls']
        for call in calls:     
            smiles = call[0]  
            try:
                mol = Chem.MolFromSmiles(smiles)
                return_value.append(str(qed.default(mol)))
            except:
                return_value.append("")

        return_json = json.dumps( { "replies" :  return_value} ), 200
        return return_json
    except Exception:
        return json.dumps( { "errorMessage": 'something unexpected in input' } ), 400