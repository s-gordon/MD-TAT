# Molecular Dynamics Trajectory Analysis Tools (MD-TAT)

## Version

Beta 0.10

## Introduction and Foreword

You are free to use it, but it may not be entirely suitable for what you are
trying to achieve. Please email feedback, bugs or suggestions to:
<a href="mailto:se2gordon@students.latrobe.edu.au">se2gordon@students.latrobe.edu.au</a>.

## Dependencies

* Python (v2.7+)
* VMD (1.9.1+). This must be incorporated into your `$PATH` variable for
    `extract_data.sh` and `vmd_analyse_and_plot.py` to work properly.

## How it all works

This project directory is designed to be used at the completion of a run as a
place to consolidate and process trajectory data generated using the workflow
MD\_workflow [here](https://github.com/s-gordon/MD_workflow/tree/pythonic). This
must be cloned **under** the top directory of MD_workflow for everything to work
properly.

Trajectory processing and analysis are divided into three main phases divided
into three independent scripts: `extract_data.sh`,
`compress_and_concatenate_trajectories.py`, and `vmd_analyse_and_plot.py`. The
order in which these scripts are run is critical.

1. `extract_data.sh`

1. `compress_and_concatenate_trajectories.py`: Additional options, such as
   trajectory sub-sampling (a.k.a. stride) can be passed as optional arguments.
   For more information try `./compress_and_concatenate_trajectories.py --help`.

1. `vmd_analyse_and_plot.py`: Options can be passed as arguments. For more
   information try `./vmd_analyse_and_plot.py --help`.

## Authorship

* Shane Gordon (<a href="mailto:se2gordon@students.latrobe.edu.au">se2gordon@students.latrobe.edu.au</a>)
