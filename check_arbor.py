import arbor
import glob
import json
from pathlib import Path


def check_morphologies(models_path, morphologies_dir='./components/morphologies'):
    models_dict = json.load(open(models_path, 'r'))
    morph_dir = Path(morphologies_dir)
    for layer, layer_dict in models_dict['locations'].items():
        print(f'LAYER {layer}')
        for model_name, models_dict in layer_dict.items():
            print(f'MODEL {model_name}')
            for model in models_dict['models']:
                if 'morphology' not in model:
                    continue
                
                swc_path = morph_dir / model['morphology']
                try:
                    arbor.load_swc_arbor(swc_path).morphology
                    passed_swc = True
                except:
                    passed_swc = False
                
                try:
                    arbor.load_swc_neuron(swc_path).morphology
                    pass_neuron = True
                except:
                    pass_neuron = False
                
                print(f'{swc_path} {passed_swc or pass_neuron}')
            print()

if __name__ == '__main__':
    check_morphologies('model_props/v1_node_models.json')
