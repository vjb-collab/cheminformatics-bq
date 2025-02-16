from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors
from rdkit.Chem import rdqueries
import json

def rdkit_fingerprint_1024_int64(request):
    try:
        return_value = []
        request_json = request.get_json()
        calls = request_json['calls']
        
        for call in calls:     
            
            smiles = call[0]  
            
            try:
                
                mol = Chem.MolFromSmiles(smiles)
                
                fp_pattern_long = Chem.PatternFingerprint(mol,fpSize=1024,tautomerFingerprints=True)
                fp_pattern_long_as_binary = DataStructs.BitVectToBinaryText(fp_pattern_long)
                fp_pattern_long_as_binary_hex = fp_pattern_long_as_binary.hex()
               
                fingerprints = {}
                fingerprints["fp_pattern_long_as_binary_hex"] = fp_pattern_long_as_binary_hex

                fingerprints["fp_as_ints"] = convert_bv_to_64bit_ints(fp_pattern_long)
                
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

def convert_bv_to_64bit_ints(fp):
    ''' use a combination of RDKit+python functionality to convert a
    fingerprint into a list of 64bit ints 
    '''
    qt = DataStructs.BitVectToBinaryText(fp)
    words = [int.from_bytes(qt[i:i+8],'big', signed=True) for i in range(0,len(qt),8)]
    return words


def rdkit_fingerprint(request):
    try:
        return_value = []
        request_json = request.get_json()
        calls = request_json['calls']
        
        for call in calls:     
            
            smiles = call[0]  
            
            try:
                
                mol = Chem.MolFromSmiles(smiles)
                
                fp_pattern_long = Chem.PatternFingerprint(mol, tautomerFingerprints=True)
                fp_pattern_long_as_binary = DataStructs.BitVectToBinaryText(fp_pattern_long)
                fp_pattern_long_as_binary_hex = fp_pattern_long_as_binary.hex()
               

                fp_morgan = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius = 2, useChirality=True, nBits = 2048)
                fp_morgan_as_binary = DataStructs.BitVectToBinaryText(fp_morgan)
                fp_morgan_as_binary_hex = fp_morgan_as_binary.hex()
                

                fingerprints = {}

                fingerprints["fp_pattern_long_as_binary_hex"] = fp_pattern_long_as_binary_hex
                fingerprints["fp_morgan_as_binary_hex"] = fp_morgan_as_binary_hex

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


def rdkit_substructure_match(request):
    try:
        return_value = []
        request_json = request.get_json()
        calls = request_json['calls']
        
        for call in calls:     
            
            fragment_smiles = call[0]
            smiles = call[1]  
                        
            try:
                mol = Chem.MolFromSmiles(smiles)
                fragment_mol = Chem.MolFromSmiles(fragment_smiles)
                return_value.append(mol.HasSubstructMatch(fragment_mol, useChirality=True))
            
            except:
                return_value.append("")

        return_json = json.dumps( { "replies" :  return_value} ), 200
        return return_json
    except Exception:
        return json.dumps( { "errorMessage": 'something unexpected in input' } ), 400

def rdkit_molecular_descriptors(request):
    try:
        return_value = []
        request_json = request.get_json()
        calls = request_json['calls']
        
        for call in calls:     
            
            smiles = call[0]  
            
            try:
                
                mol = Chem.MolFromSmiles(smiles)
                
                descriptor_list = ["ExactMolWt", "FractionCSP3", "BalabanJ", "BertzCT", 'HallKierAlpha','HeavyAtomCount','HeavyAtomMolWt', 'MaxAbsPartialCharge',  'MaxPartialCharge', 'MolLogP', 'MolMR', 'MolWt','NHOHCount','NOCount','NumAliphaticCarbocycles','NumAliphaticHeterocycles','NumAliphaticRings','NumAromaticCarbocycles','NumAromaticHeterocycles','NumAromaticRings','NumHAcceptors','NumHDonors','NumHeteroatoms','NumRadicalElectrons','NumRotatableBonds','NumSaturatedCarbocycles','NumSaturatedHeterocycles','NumSaturatedRings','NumValenceElectrons']
                
                data = {}
                
                for descriptor_name in descriptor_list:
                    data[descriptor_name] = getattr(Descriptors, descriptor_name)(mol)
                
                molecular_descriptors_as_json_string = json.dumps(data)
                
                return_value.append(molecular_descriptors_as_json_string)
            
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

def rdkit_qed(request):
    from rdkit.Chem import QED
    try:
        return_value = []
        request_json = request.get_json()
        calls = request_json['calls']
        for call in calls:     
            smiles = call[0]  
            try:
                mol = Chem.MolFromSmiles(smiles)
                return_value.append(str(QED.qed(mol)))
            except:
                return_value.append("")

        return_json = json.dumps( { "replies" :  return_value} ), 200
        return return_json
    except Exception:
        return json.dumps( { "errorMessage": 'something unexpected in input' } ), 400

