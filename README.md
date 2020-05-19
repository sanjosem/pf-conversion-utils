# pf-conversion-utils
Python toolbox to read and export PF simulation files

Append the main directory to your `PYTHONPATH`. You will be able to `import pftools` and sub-modules.
Append the main directory to your `PATH`, you will be able to execute conversion scripts which inputs are in YAML format.

## Requirements
  - python >= 3.6
  - numpy, scipy, h5py, pandas, vtk, pyyaml

## Probe conversion

  - Typical input `input_post_probes.yaml` is provided in `example` directory
  - Path are defined either as absolute or relative to the script execution
  - Script convert multiple probe files that are found in the input directory

```bash
python3 read_temporal.py input_post_probes.yaml
```

## Surface conversion

  - Typical input `input_post_surface.yaml` is provided in `example` directory
  - Path are defined either as absolute or relative to the script execution
  - Script convert single file and a single frame in vtk format.

```bash
python3 script_clean_surface.py input_post_surface.yaml
```

## Volume conversion

  - Typical input `input_post_volume.yaml` is provided in `example` directory
  - Path are defined either as absolute or relative to the script execution
  - Script convert single file and a single frame in vtk format.

```bash
python3 script_clean_volume.py input_post_volume.yaml
```
