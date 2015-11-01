# Molecular Dynamics Trajectory Analysis Tools (MD-TAT)

## Version

Beta 0.10

## Introduction and Foreword

You are free to use it, but it may not be entirely suitable for what you are
trying to achieve. Please email feedback, bugs or suggestions to:
<a href="mailto:se2gordon@students.latrobe.edu.au">se2gordon@students.latrobe.edu.au</a>.

## Dependencies

* Python (v2.7+)
* MDtraj (v1.4.2)

## How it all works

This project directory is designed to be used at the completion of a run as a
place to consolidate and process trajectory data generated using the workflow
MD\_workflow [here](https://github.com/s-gordon/MD_workflow/tree/pythonic).
You'll need to point `compress.py` in the right direction of the trajectory
files relative to where you run these scripts.

Trajectory processing and analysis are divided into three main phases divided
into two independent scripts: `compress.py`, and `analyze_plot.py`. The order in
which these scripts are run is critical.

1. `compress.py`: Additional options, such as trajectory sub-sampling (a.k.a.
   stride) can be passed as optional arguments.  For more information try
   `./compress.py --help`.

1. `analyze_plot.py`: Options can be passed as arguments. For more information
   try `./analyse_plot.py --help`.

## Authorship

* Shane Gordon (<a href="mailto:se2gordon@students.latrobe.edu.au">se2gordon@students.latrobe.edu.au</a>)
