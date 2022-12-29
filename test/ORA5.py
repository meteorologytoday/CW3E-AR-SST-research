import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-oras5',
    {
        'format': 'zip',
        'product_type': 'consolidated',
        'variable': [
            'mixed_layer_depth_0_01', 'mixed_layer_depth_0_03',
        ],
        'vertical_resolution': 'single_level',
        'month': [
            '01', '02', '03',
        ],
        'year': '1958',
    },
    'download.zip')

if False:
    c.retrieve(
        'reanalysis-oras5',
        {
            'format': 'zip',
            'product_type': 'consolidated',
            'variable': [
                'potential_temperature', 'salinity',
            ],
            'vertical_resolution': 'all_levels',
            'month': [
                '01', '02', '03',
            ],
            'year': '1958',
        },
        'download.zip')
