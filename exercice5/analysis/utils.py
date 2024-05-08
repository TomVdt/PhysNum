import subprocess
import numpy as np
from typing import Any
from itertools import product
from concurrent.futures import ProcessPoolExecutor, Future

# Matplotlib stuff
rcParams = {
    "text.usetex": True,
    "font.family": "serif",
    "axes.labelsize": 14,
    "font.size": 12,
    "legend.fontsize": 11,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "figure.figsize": (4, 3),
    "text.latex.preamble": "\n".join([
        r"\usepackage[utf8]{inputenc}",
        r"\usepackage[T1]{fontenc}",
        r"\usepackage[detect-all,locale=FR]{siunitx}",
    ]),
    'lines.markersize': 8,
    'lines.color': 'grey',
    'scatter.marker': '+',
    'errorbar.capsize': 3,
    'savefig.bbox': 'tight',
    'axes.formatter.useoffset': False,
    'axes.spines.right': False,
    'axes.spines.top': False,
}

# Simulation and paths stuff
path = '../'
executable = 'bin/ex5'
export_path = path + 'rapport/figures/'
data_path = 'data/'
config_path = 'bin/'
config_ext = '.conf'

# Utils
def stringify_dict(d: dict, sep=',') -> str:
    return sep.join(map(lambda a: str(a[0]) + "=" + str(a[1]), tuple(d.items())))

def gen_variations(kwargs: dict[str, Any]) -> list[dict[str, Any]]:
    return list(
        {a: b for a, b in zip(kwargs.keys(), c)} for c in product(*kwargs.values())
    )

# Simulations
def run(config_file: str, output_file: str, params: dict = {}) -> None:
    options = stringify_dict(params, sep=' ')
    cmd = f"{path}{executable} {path}{config_file} output='{path}{output_file}' {options}"
    # print(f"Running command `{cmd}`\n", end='')
    subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)

def run_with_params(config_name: str, all_params: list[dict[str, Any]]) -> list[tuple[dict[str, Any], np.ndarray]]:
    count = 0
    def done(fut: Future) -> None:
        nonlocal count
        count += 1
        print(f'\rRunning simulations... Done {count}/{len(all_params)}', end='')

    # Run simulations *IN PARALLEL*
    outputs = []
    with ProcessPoolExecutor(max_workers=8) as p:
        for params in all_params:
            options = stringify_dict(params)
            output_file = f"{data_path}{config_name}{',' if options else ''}{options}"
            outputs.append(output_file)
            future = p.submit(run, f'{config_path}{config_name}{config_ext}', output_file, params)
            future.add_done_callback(done)
    print()

    dataset = []
    for file, params in zip(outputs, all_params):
        data_x = np.loadtxt(path + file + "_x.out")
        data_v = np.loadtxt(path + file + "_v.out")
        data_f = np.loadtxt(path + file + "_f.out")
        data_h0 = np.loadtxt(path + file + "_h0.out")
        dataset.append((params, data_x, data_v, data_f, data_h0))
    return dataset

def load_conf(config_name: str) -> dict[str, Any]:
    conf = {}
    with open(path + config_path + config_name + config_ext, 'r') as f:
        lines = f.read().split('\n')

    for line in lines:
        if not line or line.startswith(('%', '#', '//')) or not line.strip():
            continue
        name, _, val, *_ = line.split(' ')
        name = name.strip()
        val = val.strip()
        try:
            conf[name] = float(val)
        except ValueError:
            conf[name] = val
    
    return conf
