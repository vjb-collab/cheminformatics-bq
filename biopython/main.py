import json

def biopython_sequence_complement(request):
    from Bio.Seq import Seq
    try:
        return_value = []
        request_json = request.get_json()
        calls = request_json['calls']
        
        for call in calls:     
            sequence = call[0]  
            try:
                seq = Seq(sequence)
                return_value.append( str(seq.complement()) )
            except:
                return_value.append("")

        return_json = json.dumps( { "replies" :  return_value} ), 200
        return return_json
    except Exception:
        return json.dumps( { "errorMessage": 'something unexpected in input' } ), 400

def biopython_sequence_translate(request):
    from Bio.Seq import Seq
    try:
        return_value = []
        request_json = request.get_json()
        calls = request_json['calls']
        
        for call in calls:     
            sequence = call[0]  
            try:
                seq = Seq(sequence)
                return_value.append( str(seq.translate()) )
            except:
                return_value.append("")

        return_json = json.dumps( { "replies" :  return_value} ), 200
        return return_json
    except Exception:
        return json.dumps( { "errorMessage": 'something unexpected in input' } ), 400

def biopython_sequence_translate_to_stop(request):
    from Bio.Seq import Seq
    try:
        return_value = []
        request_json = request.get_json()
        calls = request_json['calls']
        
        for call in calls:     
            sequence = call[0]  
            try:
                seq = Seq(sequence)
                return_value.append( str(seq.translate(to_stop=True)) )
            except:
                return_value.append("")

        return_json = json.dumps( { "replies" :  return_value} ), 200
        return return_json
    except Exception:
        return json.dumps( { "errorMessage": 'something unexpected in input' } ), 400