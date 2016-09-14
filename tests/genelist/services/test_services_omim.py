from genelist.services.omim import OMIM
import yaml

def init(config_stream):
    config = yaml.load(config_stream)
    api_key = config['OMIM']['api_key']
    return OMIM(api_key=api_key)

def test_parse_phenotypic_disease_models_ext(config_stream):
    omim = init(config_stream)

    gene = omim.gene(mim_number='609300')
    models = omim.parse_phenotypic_disease_models_ext(gene['phenotypes'])

    expected_models = {202110: [
        {
            'description': '17,20-lyase deficiency, isolated',
            'models': ['AR']
        },
        {
            'description': '17-alpha-hydroxylase/17,20-lyase deficiency',
            'models': ['AR']
        }
    ]}

    assert models == expected_models

    gene = omim.gene(mim_number='606609')
    models = omim.parse_phenotypic_disease_models_ext(gene['phenotypes'])

    expected_models = {
        192315: [
        {
            'description': 'Vasculopathy, retinal, with cerebral leukodystrophy',
            'models': ['AD']
        }],
        225750: [
        {
            'description': 'Aicardi-Goutieres syndrome 1, dominant and recessive',
            'models': ['AD', 'AR']
        }],
        152700: [
            {
                'description': '{Systemic lupus erythematosus, susceptibility to}',
                'models': ['AD']
        }],
        610448: [
            {
                'description': 'Chilblain lupus',
                'models': ['AD']
            }]
        }


    assert models == expected_models
