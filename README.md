# ePICRecoDataAnalysis
Analysis of ePIC simulated data. The code is written so that it can be easily converted between 3 different data reading modes:
- ROOT trees
- `podio::ROOTReader`
- `podio::ROOTFrameReader` (to be deprecated in podio v1.0)
- `podio::EventStore` (abandoned)


## Simplified version

To get a simple branch use the command below. The simple branch is a good starting point.

```Sh
git checkout simple
```


### To run interactively:


```Sh
root -l -b -q readFrameRoot.C+ | tee run.log
```

Old version:

```Sh
root -l -b -q readTreeSim.C+ | tee run.log
```

### To compile executable:

```Sh
make
```

and Run:

```Sh
./readFrameRootMain | tee run.log
```

### Batch scripts to use condor on RCF

```Sh
submitSimAnalysis.job
runSimAnaBatch.sh
```
Submit with:

```Sh
condor_submit submitSimAnalysis.job
```

### Python version

You can use `ipython3` to have interactive `python` and autocomplete feauture

```Sh
python readEvents.py
```

### TO DO
1. Add more histograms and functionality
2. Add macros for drawing histograms