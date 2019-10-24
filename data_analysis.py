import re, pickle
import matplotlib.pyplot as plt
from tabulate import tabulate
from statistics import mean, stdev
import math


def make_table():
    def avg_swap_to_string(avg_swap):
        if avg_swap is None:
            return
        return f'{avg_swap:.1f}'


    def avg_time_to_string(avg_time):
        if avg_time is None:
            return
        if avg_time >= 1:
            return f'{avg_time:.1f}'
        nzeros = -math.log10(avg_time)
        ndigits = math.floor(nzeros + 1)
        return '{0:.{1}f}'.format(avg_time, ndigits)

    spect_avg_swap, spect_stdev_swap, spect_avg_time, spect_stdev_time = get_averages_and_stdev(SPECTRAL_FILES)
    arct_greedy_avg_swap, arct_greedy_stdev_swap, arct_greedy_avg_time, arct_greedy_stdev_time = get_averages_and_stdev(ARCT_GREEDY_FILES)
    arct_extend_avg_swap, arct_extend_stdev_swap, arct_extend_avg_time, arct_extend_stdev_time = get_averages_and_stdev(ARCT_EXTEND_FILES)
    arct_simple_avg_swap, arct_simple_stdev_swap, arct_simple_avg_time, arct_simple_stdev_time = get_averages_and_stdev(ARCT_SIMPLE_FILES)
    arct_qiskit_avg_swap, arct_qiskit_stdev_swap, arct_qiskit_avg_time, arct_qiskit_stdev_time = get_averages_and_stdev(ARCT_QISKIT_FILES)
    ibm_qx_avg_swap, ibm_qx_stdev_swap, ibm_qx_avg_time, ibm_qx_stdev_time = get_averages_and_stdev(IBM_QX_MAPPING_FILES)

    values = []
    with open('benchmarks_run.txt', 'r') as f:
        sorted_files_text = f.readlines()

    for i in range(1, len(sorted_files_text)):
        line = sorted_files_text[i]
        (filename, ncnots, nqubits) = re.findall("(.*)\.qasm,(.*),(.*)", line)[0]
        new_value = [filename, nqubits, ncnots,
                     avg_swap_to_string(spect_avg_swap[filename]), avg_time_to_string(spect_avg_time[filename]),
                     avg_swap_to_string(arct_greedy_avg_swap[filename]), avg_time_to_string(arct_greedy_avg_time[filename]),
                     avg_swap_to_string(arct_simple_avg_swap[filename]), avg_time_to_string(arct_simple_avg_time[filename]),
                     avg_swap_to_string(arct_extend_avg_swap[filename]), avg_time_to_string(arct_extend_avg_time[filename]),
                     avg_swap_to_string(arct_qiskit_avg_swap[filename]), avg_time_to_string(arct_qiskit_avg_time[filename]),
                     avg_swap_to_string(ibm_qx_avg_swap.get(filename, None)), avg_time_to_string(ibm_qx_avg_time.get(filename, None))]
        values.append(new_value)
    headers = ['\\textsc{Benchmark Circuit}',
               '\\textsc{Qubits}',
               '\\textsc{CNOTs}',
               '\\textsc{SWAP Count}',
               '\\textsc{CPU Time (sec)}',
               '\\textsc{SWAP Count}',
               '\\textsc{CPU Time (sec)}',
               '\\textsc{SWAP Count}',
               '\\textsc{CPU Time (sec)}',
               '\\textsc{SWAP Count}',
               '\\textsc{CPU Time (sec)}',
               '\\textsc{SWAP Count}',
               '\\textsc{CPU Time (sec)}',
               '\\textsc{SWAP Count}',
               '\\textsc{CPU Time (sec)}']

    return tabulate(values, headers=headers, tablefmt="latex", missingval="out of memory")


def get_averages_and_stdev(summary_files):
    directory = 'results/summary/'
    swap_data = {}
    time_data = {}

    for summary_file in summary_files:
        with open(directory + summary_file, 'r') as f:
            lines = f.readlines()

        for line in lines:
            regex = "(.*)\.qasm,(.*),(.*)"
            res = re.findall(regex, line)
            if not res:
                continue
            filename, swap_count, time = res[0]
            swap_data.setdefault(filename, []).append(int(swap_count))
            time_data.setdefault(filename, []).append(float(time))

    for key in set(swap_data.keys()):
        if len(swap_data[key]) < 5:
            del swap_data[key]
            del time_data[key]
    avg_swap = {key: mean(swap_data[key]) for key in swap_data}
    stdev_swap = {key: stdev(swap_data[key]) for key in swap_data}
    avg_time = {key: mean(time_data[key]) for key in time_data}
    stdev_time = {key: stdev(time_data[key]) for key in time_data}

    return avg_swap, stdev_swap, avg_time, stdev_time


def ratio_data(summaries1, summaries2):
    avg_swap_1, stdev_swap_1, avg_time_1, stdev_time_1 = get_averages_and_stdev(summaries1)
    avg_swap_2, stdev_swap_2, avg_time_2, stdev_time_2 = get_averages_and_stdev(summaries2)
    swap_ratios = {}
    swap_errors = {}
    time_ratios = {}
    time_errors = {}
    testfiles = set(avg_swap_1.keys()) | set(avg_swap_2.keys())
    for testfile in testfiles:
        if testfile not in avg_swap_1 or testfile not in avg_swap_2:
            continue
        swap_ratios[testfile] = avg_swap_1[testfile] / avg_swap_2[testfile] \
            if avg_swap_1[testfile] + avg_swap_2[testfile] > 0 else 1 # hack for both being 0
        time_ratios[testfile] = avg_time_1[testfile] / avg_time_2[testfile]
        swap_errors[testfile] = 0 if avg_swap_1[testfile] == 0 \
            else swap_ratios[testfile] * ((stdev_swap_1[testfile] / avg_swap_1[testfile])**2
                                 + (stdev_swap_2[testfile] / avg_swap_2[testfile])**2)**0.5
        time_errors[testfile] = time_ratios[testfile] * ((stdev_time_1[testfile] / avg_time_1[testfile])**2
                                 + (stdev_time_2[testfile] / avg_time_2[testfile])**2)**0.5
    return swap_ratios, swap_errors, time_ratios, time_errors


def make_summary_ratio_plots(summaries1, summaries2, label1, label2,
                             annotations, titlelabel2=None, xscalelog=False):
    swap_ratios, swap_errors, time_ratios, time_errors = ratio_data(summaries1, summaries2)
    x = []
    y = []
    xerr = []
    yerr = []
    
    x_2 = []
    y_2 = []
    xerr_2 = []
    yerr_2 = []
    for filename in swap_ratios:
        y.append(swap_ratios[filename])
        yerr.append(swap_errors[filename])
        x.append(time_ratios[filename])
        xerr.append(time_errors[filename])
        if 'after' in filename:
            y_2.append(swap_ratios[filename])
            yerr_2.append(swap_errors[filename])
            x_2.append(time_ratios[filename])
            xerr_2.append(time_errors[filename])

    if xscalelog:
        plt.xscale('log')

    plt.errorbar(x, y,
                 xerr=xerr,
                 yerr=yerr,
                 fmt='o',
                 markersize=10,
                 elinewidth=3,
                 capsize=3,
                 zorder=2)

    plt.hlines([1], -0.05, max(pt + err for pt, err in zip(x, xerr))*1.1, linewidth=5, linestyle='--', zorder=1)
    plt.vlines([1], -0.05, max(pt + err for pt, err in zip(y, yerr)), linewidth=5, linestyle='--', zorder=1)
    plt.xlabel(f"CPU Time Ratios ({label1} / {label2})", fontsize=25)
    plt.ylabel(f"Swap Count Ratios ({label1} / {label2})", fontsize=25)
    plt.title(f"Swap Count and Runtime Ratios of\n{label1} / {titlelabel2 if titlelabel2 else label2 }",
              fontsize=25)
    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.xlim(xmax=max(pt + err for pt, err in zip(x, xerr))*1.2)

    for filename in annotations:
        plt.annotate(filename, (time_ratios[filename], swap_ratios[filename]),
                     xytext=annotations[filename],
                     arrowprops=dict(arrowstyle='->', color='red'),
                     fontsize=25)
    plt.show()


SPECTRAL_FILES = ['timing_greedy_spectral_summary_A2_v2_2019_07_09.txt',
                  'timing_greedy_spectral_summary_A2_v2_2019_07_10.txt',
                  'timing_greedy_spectral_summary_A2_v2_2019_07_10_trial_2.txt',
                  'timing_greedy_spectral_summary_A2_v2_2019_07_10_trial_3.txt',
                  'timing_greedy_spectral_summary_A2_v2_2019_07_19.txt',
                  'timing_greedy_spectral_summary_A2_v2_2019_08_05.txt',
                  'timing_greedy_spectral_summary_A2_v2_2019_08_06.txt',
                  'timing_greedy_spectral_summary_A2_v2_2019_08_06_trial_2.txt',
                  'timing_greedy_spectral_summary_A2_v2_2019_08_06_trial_3.txt',
                  'timing_greedy_spectral_summary_A2_v2_2019_08_06_trial_4.txt',]

ARCT_GREEDY_FILES = ['timing_arct_summary_A2_v2_greedy.txt',
                     'timing_arct_summary_A2_v2_greedy_2019_07_10.txt',
                     'timing_arct_summary_A2_v2_greedy_2019_07_12.txt',
                     'timing_arct_summary_A2_v2_greedy_2019_07_17.txt',
                     'timing_arct_summary_A2_v2_greedy_2019_07_20.txt',
                     'timing_arct_summary_A2_v2_greedy_2019_08_06_trial_1.txt',
                     'timing_arct_summary_A2_v2_greedy_2019_08_06_trial_2.txt',
                     'timing_arct_summary_A2_v2_greedy_2019_08_06_trial_3.txt',
                     'timing_arct_summary_A2_v2_greedy_2019_08_06_trial_4.txt',
                     'timing_arct_summary_A2_v2_greedy_2019_08_06_trial_5.txt', ]

ARCT_EXTEND_FILES = ['timing_arct_summary_A2_v2_extend.txt',
                     'timing_arct_summary_A2_v2_extend_2019_07_10.txt',
                     'timing_arct_summary_A2_v2_extend_2019_07_12.txt',
                     'timing_arct_summary_A2_v2_extend_2019_07_17.txt',
                     'timing_arct_summary_A2_v2_extend_2019_07_20.txt',
                     'timing_arct_summary_A2_v2_extend_2019_08_06_trial_1.txt',
                     'timing_arct_summary_A2_v2_extend_2019_08_06_trial_2.txt',
                     'timing_arct_summary_A2_v2_extend_2019_08_06_trial_3.txt',
                     'timing_arct_summary_A2_v2_extend_2019_08_06_trial_4.txt',
                     'timing_arct_summary_A2_v2_extend_2019_08_06_trial_5.txt',]

ARCT_SIMPLE_FILES = ['timing_arct_summary_A2_v2_simple.txt',
                     'timing_arct_summary_A2_v2_simple_2019_07_15.txt',
                     'timing_arct_summary_A2_v2_simple_2019_07_16.txt',
                     'timing_arct_summary_A2_v2_simple_2019_07_17.txt',
                     'timing_arct_summary_A2_v2_simple_2019_07_22.txt',
                     'timing_arct_summary_A2_v2_simple_2019_08_06_trial_1.txt',
                     'timing_arct_summary_A2_v2_simple_2019_08_06_trial_2.txt',
                     'timing_arct_summary_A2_v2_simple_2019_08_06_trial_3.txt',
                     'timing_arct_summary_A2_v2_simple_2019_08_06_trial_4.txt',
                     'timing_arct_summary_A2_v2_simple_2019_08_06_trial_5.txt', ]

ARCT_QISKIT_FILES = ['timing_arct_summary_A2_v2_qiskit.txt',
                     'timing_arct_summary_A2_v2_qiskit_2019_07_15.txt',
                     'timing_arct_summary_A2_v2_qiskit_2019_07_16.txt',
                     'timing_arct_summary_A2_v2_qiskit_2019_07_17.txt',
                     'timing_arct_summary_A2_v2_qiskit_2019_07_22.txt',
                     'timing_arct_summary_A2_v2_qiskit_2019_08_06_trial_1.txt',
                     'timing_arct_summary_A2_v2_qiskit_2019_08_06_trial_2.txt',
                     'timing_arct_summary_A2_v2_qiskit_2019_08_06_trial_3.txt',
                     'timing_arct_summary_A2_v2_qiskit_2019_08_06_trial_4.txt',
                     'timing_arct_summary_A2_v2_qiskit_2019_08_06_trial_5.txt', ]

IBM_QX_MAPPING_FILES = ['timing_ibm_qx_mapping_summary_A2_v2.txt',
                        'timing_ibm_qx_mapping_summary_A2_v2_2019_07_15_trial_1.txt',
                        'timing_ibm_qx_mapping_summary_A2_v2_2019_07_15_trial_2.txt',
                        'timing_ibm_qx_mapping_summary_A2_v2_2019_07_23_trial_1.txt',
                        'timing_ibm_qx_mapping_summary_A2_v2_2019_07_23_trial_2.txt',
                        'timing_ibm_qx_mapping_summary_A2_v2_2019_08_06_trial_1.txt',
                        'timing_ibm_qx_mapping_summary_A2_v2_2019_08_06_trial_2.txt',
                        'timing_ibm_qx_mapping_summary_A2_v2_2019_08_06_trial_3.txt',
                        'timing_ibm_qx_mapping_summary_A2_v2_2019_08_06_trial_4.txt',
                        'timing_ibm_qx_mapping_summary_A2_v2_2019_08_06_trial_5.txt',]

if __name__ == '__main__':
    with open('table.tex', 'w') as f:
        f.write(make_table())

    make_table()

    spectral_label = "Spectral"
    arct_label = "CSU19"
    arct_greedy_title = "CSU19 Greedy Size"
    arct_extend_title = "CSU19 Extend Size"
    arct_simple_title = "CSU19 Simple Size"
    arct_qiskit_title = "CSU19 Qiskit-Based"
    ibm_qx_label = "ZPW18"

    ibm_qx_annotations = {'graycode6_47': (40, 0.2),
                          'tof_10_after_heavy': (40, 4.5),
                          'tof_10_after_light': (200, 4)}

    make_summary_ratio_plots(SPECTRAL_FILES, IBM_QX_MAPPING_FILES,
                             spectral_label, ibm_qx_label,
                             ibm_qx_annotations, xscalelog=True)

    arct_greedy_annotations = {"ising_model_10": (1.5, 0.2),
                               "ising_model_13": (1.5, 0.1),
                               "ising_model_16": (1.5, 0),
                               "graycode6_47": (1.5, 0.3),
                               "life_238": (2.50, 0.6)}
    make_summary_ratio_plots(SPECTRAL_FILES, ARCT_GREEDY_FILES,
                             spectral_label, arct_label,
                             arct_greedy_annotations, titlelabel2=arct_greedy_title)

    arct_extend_annotations = {"ising_model_10": (0.3, 0.3),
                               "ising_model_13": (0.1, 0.2),
                               "ising_model_16": (0.1, 0.1),
                               "graycode6_47": (0.3, 0.0),
                               "barenco_tof_4_after_heavy": (1.4, 1.1),
                               "barenco_tof_4_after_light": (1.5, 0.6),
                               }
    make_summary_ratio_plots(SPECTRAL_FILES, ARCT_EXTEND_FILES,
                             spectral_label, arct_label,
                             arct_extend_annotations, titlelabel2=arct_extend_title)

    arct_simple_annotations = {"rd84_253": (1.5, 1.5),
                               "ising_model_10": (0.4, 0.36),
                               "ising_model_13": (0.4, 0.24),
                               "ising_model_16": (0.4, 0.12),
                               "graycode6_47": (0.4, 0.0),
                               "barenco_tof_4_after_heavy": (1.6, 1.2),
                               "barenco_tof_4_after_light": (1.8, 0.7),}
    make_summary_ratio_plots(SPECTRAL_FILES, ARCT_SIMPLE_FILES,
                             spectral_label, arct_label,
                             arct_simple_annotations, titlelabel2=arct_simple_title)

    arct_qiskit_annotations = {"ising_model_10": (1.4, 0.43),
                               "ising_model_13": (1.4, 0.29),
                               "ising_model_16": (1.4, 0.15),
                               "graycode6_47": (1.4, 0.0),
                               "life_238": (2.5, 1.1),
                               "sym9_148": (2.5, 0.5)}
    make_summary_ratio_plots(SPECTRAL_FILES, ARCT_QISKIT_FILES,
                             spectral_label, arct_label,
                             arct_qiskit_annotations, titlelabel2=arct_qiskit_title)
