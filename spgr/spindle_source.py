
hemi_subj = {'EM09': 'rh',
             'MG17': 'rh',
             'MG33': 'lh',
             'MG37': 'lh',
             'MG61': 'lh',
             'MG63': 'rh',
             }



def get_morph_linear(subj, data_options):

    chan = get_chan_used_in_analysis(subj, 'sleep', ('grid', ),
                                     **data_options)[1]
    fs = Freesurfer(join(REC_DIR, subj, FS_PATH))
    surf = fs.read_surf(hemi_subj[subj])
    l = Linear(surf, chan, std=STD, threshold=THRESHOLD)

    if hemi_subj[subj] == 'lh':
        for one_chan in chan.chan:
            one_chan.xyz *= (-1, 1, 1)
    surf = fs.read_surf('rh')
    l = Linear(surf, chan, std=STD, threshold=THRESHOLD)
    m = Morph(surf)

    return