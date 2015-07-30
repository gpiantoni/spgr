from collections import OrderedDict

from numpy import arange, zeros

from .constants import avg_regions

LOBES = (('frontal', ['caudalanteriorcingulate',
                      'caudalmiddlefrontal',
                      'lateralorbitofrontal',
                      'medialorbitofrontal',
                      'paracentral',
                      'parsopercularis',
                      'parsorbitalis',
                      'parstriangularis',
                      'precentral',
                      'rostralanteriorcingulate',
                      'rostralmiddlefrontal',
                      'superiorfrontal',
                      'frontalpole',
                      ]),
         ('parietal', ['inferiorparietal',
                       'isthmuscingulate',
                       'postcentral',
                       'posteriorcingulate',
                       'precuneus',
                       'superiorparietal',
                       'supramarginal',
                       ]),
         ('temporal', ['bankssts',
                       'entorhinal',
                       'fusiform',
                       'inferiortemporal',
                       'middletemporal',
                       'parahippocampal',
                       'superiortemporal',
                       'transversetemporal',
                       'temporalpole',
                       ]),
         ('occipital', ['cuneus',
                        'lateraloccipital',
                        'lingual',
                        'pericalcarine',
                        ]),
         ('insula', ['insula',
                     ]),
         )
LOBES = OrderedDict(LOBES)

LOBE_COLOR_LIMITS = (0, 66)

LOBE_COLORS = {'frontal': 60,
               'parietal': 24,
               'occipital': 41,
               'temporal': 10,
               'insula': 32,
               'unknown': 32}


def freesurfer_color_code():
    TOTAL_N_REGIONS = 36
    freesurfer_code = arange(TOTAL_N_REGIONS)
    freesurfer_code[0] = -1

    colorcode = zeros(TOTAL_N_REGIONS)
    for i, region in enumerate(avg_regions):
        lobe = _find_lobe(region)
        colorcode[i] = LOBE_COLORS[lobe]

    return freesurfer_code, colorcode


def _find_lobe(region):
    for lobe_name, lobe_areas in LOBES.items():
        if region in lobe_areas:
            return lobe_name

    return 'unknown'
