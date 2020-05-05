import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.io import loadmat, savemat
import pickle
import sys
import argparse
import os

try:
    path2spikefinder = 'spikefinder-python-master/'
    sys.path.append(path2spikefinder)
    import spikefinder
except:
    print('Need to download SpikeFinder from https://github.com/codeneuro/spikefinder-python')



kind_interpolation = 'next'


def load_data_set(path, folder, dataset, train=True):
    if train:
        calcium = pd.read_csv(path + folder + dataset +
                              '.train.calcium.nopreprocessing.csv')
        spikes = pd.read_csv(path + folder + dataset +
                             '.train.spike.nopreprocessing.csv')
    else:
        calcium = pd.read_csv(path + folder + dataset +
                              '.test.calcium.nopreprocessing.csv')

    all_ids = []
    all_T = []
    all_dt = []
    all_F = []
    if train:
        all_N = []

    for index in calcium.columns:
        all_ids.append(index)
        F = np.array(calcium[index])
        dt = 1.0 / F[0]
        F = F[1:]
        F = F[~np.isnan(F)]
        all_T.append(len(F))
        all_F.append(F.copy())
        all_dt.append(dt)

        if train:
            Nraw = np.array(spikes[index])
            Nraw = Nraw[~np.isnan(Nraw)]
            N = np.zeros(F.shape)
            for t in Nraw:
                time_discrete = int(np.ceil((t / 1000) / dt)) + 1 - 1
                # i) mesure en ms
                # ii) je rajoute 1 car la premi??re mesure a lieu ?? t=0 au lieu de t=dt.
                # iii) j'enleve -1 car les indices commencent ?? z??ro.
                # Exemple: t/1000 = 0.1 * dt. => Le spike a eu lieu entre la premi??re et la deuxi??me mesure => spike[1] = 1.
                if (time_discrete >= 0) & (time_discrete < len(F)):
                    N[time_discrete] = 1
            all_N.append(N)

    if min(all_T) == max(all_T):  # Otherwise, all_F is
        all_F.append([])
        all_N.append([])

    env = {}
    env['ids'] = np.array(all_ids, dtype=np.object)
    env['dt'] = np.array(all_dt)
    env['T'] = np.array(all_T)
    env['F'] = np.array(all_F, dtype=np.object)
    env['nNeurons'] = len(calcium.columns)
    if train:
        env['N'] = np.array(all_N, dtype=np.object)
    return env


def make_matlab_files():
    path =''
    folder = 'datasets 1-5 nopreprocessing/'

    for dataset in ['1', '2', '3', '4', '5']:
        for train in [True, False]:
            env = load_data_set(path, folder, dataset, train=train)
            if train:
                savemat(path + 'matlab_processed/%s.train.mat' % dataset, env)
            else:
                savemat(path + 'matlab_processed/%s.test.mat' % dataset, env)
    return


def load_matlab_results(folder, dataset, train, tauBaseline=None):
    if train & (tauBaseline is not None):
        name = '%s.train_baseline_%s.mat' % (dataset, tauBaseline)
    elif (not train) & (tauBaseline is not None):
        name = '%s.test_baseline_%s.mat' % (dataset, tauBaseline)
    if train & (tauBaseline is None):
        name = '%s.train.mat' % dataset
    elif (not train) & (tauBaseline is None):
        name = '%s.test.mat' % dataset
    env = loadmat(path + folder + name)
    return env


def get_target_size(dataset, train):
    if train:
        folder = 'spikefinder.train/'
        prefix = '.train.'
    else:
        folder = 'spikefinder.test/'
        prefix = '.test.'
    # path = '/Volumes/Carte_SD/SpikeFinder contest/'
    path = ''
    calcium = pd.read_csv(path + folder + dataset + prefix + 'calcium.csv')
    return len(calcium[calcium.columns[0]])


def resample_predictions(N, time_matlab, time_target):
    N[0] *=0 # Remove first spike: it represents the cumulative calcium transient before the start of the recording.
    Ntarget = interp1d(time_matlab, N, fill_value=np.nan,
                       bounds_error=False, kind=kind_interpolation)(time_target)
    return Ntarget


def format_matlab_results(folder, dataset, train, offset=0, tauBaseline=None, write=True, SR=1, convolve=False):
    try:
        env = load_matlab_results(
            folder, dataset, train, tauBaseline=tauBaseline)
    except:
        print('Matlab file does not exist:',
              folder, dataset, train, tauBaseline)
        return -1

    nNeurons = int(env['nNeurons'])
    dt = env['dt']
    if dataset in ['7', '8']:
        dt = np.ones([1, nNeurons]) * 0.01665


    T_target = get_target_size(dataset, train)
    dt_target = 0.01
    time_target = np.arange(T_target) * dt_target

    Ninf_formatted = np.zeros([nNeurons, T_target])
    for i in range(nNeurons):
        dt_ = dt[0,i]
        T_ = len( env['F'][0,i][0] )        
        Ninf_matlab = np.squeeze(env['Ninf'][0,i])

        if dataset in ['6','7','8','9','10']:
            if SR == 1:
                time_matlab = env['Time'][0, i][0] - env['Time'][0, i][0][0]+ offset * dt_target
            else:
                time_matlab = interp1d( np.arange(1, T_+1) , env['Time'][0, i][0] - env['Time'][0, i][0][0], bounds_error=False,fill_value='extrapolate' )( np.arange(1, T_ * SR+1)/SR ) + offset * dt_target
        else:
            if SR == 1:
                time_matlab = np.arange(T_) * dt_  + offset * dt_target
            else:
                time_matlab = np.arange(T_) * dt_/SR - 0.5 * dt_/SR  + offset * dt_target



        Ninf_formatted[i] = resample_predictions(Ninf_matlab, time_matlab, time_target)
        if convolve:
            PSF = env['all_PSF'][i, :]
            PSF[np.isnan(PSF)] = 0
            Ninf_formatted[i] = np.convolve(
                PSF, Ninf_formatted[i], 'same')
    DF_Ninf_formatted = pd.DataFrame(data=Ninf_formatted.T, columns=[
                                     '%s' % i for i in range(nNeurons)])

    if write:
        try:
            os.mkdir(path + folder + '/csv_final5/')
        except:
            pass
        try:
            os.mkdir(path + folder + '/csv_final_convolve5/')
        except:
            pass
        if train & convolve:
            name = path + folder + \
                '/csv_final_convolve5/%s.train.spikes.csv' % (dataset)
        elif train & (not convolve):
            name = path + folder + \
                '/csv_final5/%s.train.spikes.csv' % (dataset)
        elif (not train) & convolve:
            name = path + folder + \
                '/csv_final_convolve5/%s.test.spikes.csv' % (dataset)
        elif (not train) & (not convolve):
            name = path + folder + \
                '/csv_final5/%s.test.spikes.csv' % (dataset)

        DF_Ninf_formatted.to_csv(name, sep=',',float_format='%.4f',index=False)
    return DF_Ninf_formatted


# Je veux une fonction qui:
# Pour chaque methode (= pour un folder donn??):
# Pour chaque dataset:
# Charge l'offset par neurone bas?? sur le d??calage entre mes spikes ground truth et ceux de la comp??tition.
# Pour chaque baseline:
# Pour chaque offset entre -8 dt et 8 dt (li?? au d??calage entre fluorescence et spikes ground truth):
# Charge les r??sutlats matlab, les formatte ?? 100 Hz et r??ecrit en csv.
# Charge le csv et fait la comparaison de scoring.
# Pour chaque dataset, calcule la combinaison d'offset et de baseline optimale.
# Pour chaque dataset, copie le fichier correspondant dans spike.train. ....
# Retourne le score, l'offset, et la baseline pour chaque dataset.

def analyze_results(folder, datasets=[str(x) for x in range(1, 11)], train_test=[True, False], offsets=np.zeros(10),
                    tauBaselines=[10, 10,10,10,10,10], SRs=[1,1,1,1,1,1,1,1,1,1], convolve=False):
    if type(datasets) is not list:
        datasets = [datasets]
    if type(train_test) is not list:
        train_test = [train_test]

    if offsets is None:
        if SR > 1:
            offsets = np.arange(-10, 10 + 0.5, 0.5)
        else:
            offsets = np.arange(-10, 10 + 1.0, 1.0)

    ndatasets = len(datasets)

    correlations = np.zeros([ndatasets])

    for i, dataset in enumerate(datasets):
        for train in train_test:
            tauBaseline = tauBaselines[i]
            offset = offsets[i]
            print('dataset: %s, train: %s, tauBaseline: %s s, offset: %s' % (
                dataset, train, tauBaseline, offset))
            Ninf_formatted = format_matlab_results(folder,
                                                   dataset, train, offset=offset, tauBaseline=tauBaseline, SR=SRs[i], convolve=convolve)
            if train:
                Ntarget = pd.read_csv(
                    path + 'spikefinder.train/' + dataset + '.train.spikes.csv')                    
                correlations[i] = np.median(spikefinder.score(
                    Ntarget, Ninf_formatted, method='corr'))

    return correlations


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', type=str, default='matlab_results/',
                        help='Which folder to analyze')
    parser.add_argument('--SR', type=int, default=0,
                        help='Superresolution?')
    parser.add_argument('--convolve', type=bool, default=0,
                        help='Convolve with PSF?')

    FLAGS = parser.parse_args()
    if FLAGS.SR == 1:
        SRs = [1,1,1,1,1,2,2,2,6,3]
    else:
        SRs = [1,1,1,1,1,1,1,1,1,1]

    tauBaselines = [60, 10, 60, 30, 40, 40, 10, 60, 30, 20]
    # offsets = [0,-2,2,1,-9, 0,0,0,0,0]
    offsets = [0,0,0,0,0,0,0,0,0,0]
    datasets = [str(x) for x in range(1, 11)]

    if FLAGS.SR == 1:
        tauBaselines = tauBaselines[5:]
        offsets = offsets[5:]
        datasets = datasets[5:]
        SRs = SRs[5:]


    correlations= analyze_results(
        FLAGS.folder, SRs= SRs,
         offsets=offsets, convolve=False,
        tauBaselines = tauBaselines, train_test = True,
        datasets=datasets)
    print('correlations')
    print(correlations)
    env = {'tauBaselines':tauBaselines,'offsets':offsets,'datasets':datasets,'correlations': correlations}

    pickle.dump(env, open(path + FLAGS.folder + 'analyzed_results5_%s_%s.data' %
                          (FLAGS.convolve, kind_interpolation), 'wb'))

    print('Done!')
